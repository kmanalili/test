# Standard library
import os
import glob
import time
import re
import warnings
import subprocess
import itertools

# Scientific and plotting
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from skimage.draw import line_nd
import scipy
from scipy.spatial import cKDTree
from operator import itemgetter

# Visualization
import pyvista as pv
from stl import mesh
import pymeshfix as mf

# Domain packages
import srtm
import wellpathpy as wp
from pyCRGI.pure import get_value

# Windows
try:
    import pyautogui
    import win32gui
except ImportError:
    print('Ubuntu import error, use pyautogui through Windows')

class SonarPy:

    """
    A Python Class for processing and analyzing sonar data.

    ## Parameters ---------------------------------------------------------------------------------

    surveyxyz : tuple or None
        The (x, y, z) of the survey datum in the input CRS.
    surveyxyzWGS84 : tuple or None
        The (lat, long, elevation) of the survey datum in WGS84.
    year : int, default=2024
        Year of the survey.
    metric : bool, default=False
        Whether elevation and distance units are in meters (True) or feet (False).
    crs : str or pyproj.crs.crs.CRS
        Cordinate Reference System (CRS) of the input data.
        Can be EPSG code string or pyproj CRS object.
    utm_shp_path : str or None optional
        Path to UTM Zone shapefile.
        If not provided, defaults to "C:/GIS/World_UTM_Grid.zip"`, or prompts the user.

    ## Attributes ---------------------------------------------------------------------------------

    crs : str or pyproj.CRS
        Output coordinate reference system.
    metric : bool
        Whether elevation and distance units are in meters (True) or feet (False).
    utm_shp_path : str
        File path to the UTM zone shapefile.
    utmzones : geopandas.GeoDataFrame
        UTM zone polygons.
    zonecrs : str
        EPSG code string representing the UTM zone of the survey point.
    magdec : float
        Magnetic declination in degrees.
    surveyxyz : tuple or None
        The (x, y, z) of the survey datum in the input CRS.
    pntgdf : geopandas.GeoDataFrame
        The survey point as a GeoDataFrame in the input CRS.
    pntgdfUTM : geopandas.GeoDataFrame
        The survey point reprojected to the UTM CRS.
    surveyxyzUTM : tuple
        The survey datum transformed into UTM coordinates.
    xyz_gdf : geopandas.GeoDataFrame
        3D point cloud in UTM coordinates.
    xyz_gdfCRS : geopandas.GeoDataFrame
        3D point cloud reprojected into the user-specified CRS.

    ## Setup --------------------------------------------------------------------------------------

    Download the IGRF13 coefficient file from:
        https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt

    Place it in:
        C:/ProgramData/anaconda3/envs/<yourenv>/lib/site-packages/pyCRGI/data/igrf13coeffs.txt

    ## Example ------------------------------------------------------------------------------------

    file = '2m-2204.lwt'
    surveyxyz = (-49,35,5280)
    sp_2m_2024 = SonarPy(file, surveyxyz, crs='EPSG:4326')
    xyz_gdf, xyz_gdfCRS = sp_2m_2024.generate_gdf()

    ## Goals --------------------------------------------------------------------------------------

    1) K.I.S.S.
    2) Minimize tech debt
    3) Ease of use first; high end function 2nd

    ## To Do --------------------------------------------------------------------------------------

    Add:
    -2D hull geometry output
    -automatic format detection; american, german, french formats
    -automatic assignment of the metric
    -add/transfer parser for french and german formats
    -Date is in the LWT exports from CWR files ! add

    Fix
    -round outputs in a pseudo signaficant digits
    -after using this on different file types I believe the direction of the package would be to move it toward
        standalone functions for the time being rather than an all alone process. Too much variation inbetween
        data types. In the future after there may be a shift back after the transformed longwall table

    """

    def __init__(self, surveyxyz=None, surveyxyzWGS84=None, year=2024, metric=False, crs='EPSG:4326',
                 utm_shp_path=None):
        self.__version__ = '0.0.0d'
        self.crs = crs
        self.metric = metric

        # Determine UTM shapefile path
        if utm_shp_path:
            self.utm_shp_path = utm_shp_path
        elif os.path.exists(r"C:/GIS/World_UTM_Grid.zip"):
            self.utm_shp_path = r"C:/GIS/World_UTM_Grid.zip"
        else:
            self.utm_shp_path = input('Path of World_UTM_Grid.zip?')

        # Validate that UTM shapefile exists
        if not os.path.exists(self.utm_shp_path):
            raise Exception('UTM Zone Shapefile not found in', self.utm_shp_path,
                            'Place file in C:\GIS\ or define with utm_shp_path parameter')

        # Load UTM zones shapefile
        self.utmzones = gpd.read_file(self.utm_shp_path)

        # Compute magnetic declination if WGS84 coordinates are provided
        if surveyxyzWGS84:
            d, i, h, x, y, z, f = get_value(surveyxyzWGS84[1], surveyxyzWGS84[0], surveyxyzWGS84[2] / 3280.84, year)
            self.magdec = d
        else:
            warnings.warn("Provide a surveyxyz and year to add magnetic declination correction to the sonar.")
            self.magdec = 0

        # If CRS survey point is provided, store and assign UTM zone
        if surveyxyz:
            self.surveyxyz = surveyxyz
            self.utm_epsg_code_at_point(self.surveyxyz, crs=self.crs)

        # Future extension: load sonar file
        # if self.file:
        #    self.open_lwt()
        #    self.lwt_df_to_delta_points()

    def las2txt_path(path):
        """
        Converts all LAZ files in a directory to ASCII text files using las2txt64.

        Each output file contains:
        - x, y, z coordinates (columns 1-3)
        - intensity (column 4)
        - gps_time (column 5)
        - a header row for use with txt2las64
        
        Parameters
        ----------
        path : str
            Path to the folder containing .laz files.
        """
        # Save current working directory and change to target directory
        org = os.getcwd()
        os.chdir(path)

        # Convert all .laz files to .txt using las2txt64
        os.system("las2txt64 -i *.laz -parse xyzit -sep comma")

        # Return to original working directory
        os.chdir(org)

    def add_headers_to_LAStxt(path):
        """
        Add column headers to las2txt output and save as csv.

        The input .txt files are expected to contain: X, Y, Z, intensity, and gps_time.

        Parameters
        ----------
        path : str
            Path to the folder containing .txt files exported by las2txt64.
        """
        # Ensure path ends with '/'
        path = path.replace('\\', '/')
        if path[-1] != '/':
            path = path + '/'

        # Find all txt files in directory
        files = [x.replace('\\', '/') for x in glob.glob(path + '*.txt')]

        for file in files:
            try:
                df = pd.read_csv(file, header=None)
                df.columns = ['X', 'Y', 'Z', 'intensity', 'gps_time']
                filename = file.replace('.txt', '.csv')
                df.to_csv(filename, index=False)
            except Exception as e:
                warnings.warn('ERR ' + file + ' ' + str(e))

    def dmv_exported_lwt_to_dfT(self, path):
        """
        Converts LWT files exported by DimCav Viewer to a dfT format DataFrame

        Parameters
        ----------
        path : str
            Path to the LWT file.

        Returns
        -------
        dfT : pandas.DataFrame
            Reformatted long-form DataFrame where each row represents a single sonar shot,
            including attributes such as depth, range, azimuth (cAzi), and tilt.
        """
        with open(path, 'r') as f:
            txt = f.read()

        dfT = pd.DataFrame()

        # Split into tables
        if '\n\n\n' in txt:
            tbls = txt.split('\n\n\n')
        else:
            tbls = txt.split('\n\n')

        # Parse tables into lists
        tbls = [x.split('\n') for x in tbls]
        tbls = [[y.split(';') for y in x] for x in tbls]
        tbls = [[[z for z in y if len(z) > 0] for y in x] for x in tbls]
        tbls = [x for x in tbls if len(x) > 0]
        tbls = [x for x in tbls if x != [[]]]

        for tbl in tbls:
            # Get metadata
            _, depth, _, tilt, _, vos = tbl[0]
            # Get data
            tbl = pd.DataFrame(tbl[2:], columns=tbl[1])
            tbl = tbl.astype(float)

            for index, row in tbl.iterrows():
                for col in list(tbl):
                    if col == 'Bearing':
                        bearing = row[col]
                    else:
                        plus = float(col.replace('+', '').strip())

                        r = pd.DataFrame({'depth': [depth],
                                          'tilt': [tilt],
                                          'vos': [vos],
                                          'cAzi': [bearing + plus],
                                          'r': row[col]})
                        dfT = pd.concat([dfT, r], ignore_index=True)

        dfT = dfT.astype(float)
        return dfT

    def get_well_deviation_delta(self, df, wb):
        """
        Returns the (dx, dy, dz) deviation from a wellbore survey (wb) at the correct measured depth or last depth.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing the sonar data.
        wb : pandas.DataFrame
            Wellbore survey data containing columns 'MDepth', 'dx_ft', 'dy_ft', 'dz_ft'.

        Returns
        -------
        tuple
            Delta (dx, dy, dz) in feet.
        """
        depthmin = int(df.depth.min())
        deeperI = wb[wb.MDepth >= depthmin].index.tolist()

        if len(deeperI) > 0:
            if depthmin not in wb.MDepth.values:
                depths = [x for x in range(wb.MDepth.min(), wb.MDepth.max() + 1) if x not in wb.MDepth.values]
                if depthmin not in wb.MDepth.values:
                    wb = pd.concat([wb, pd.DataFrame({'MDepth': depths})], ignore_index=True)
                    wb.sort_values(by='MDepth', inplace=True)
                wb = wb.interpolate(method='linear', axis=0)
            dx = wb[wb.MDepth == depthmin]['dx_ft'].values[0]
            dy = wb[wb.MDepth == depthmin]['dy_ft'].values[0]
            dz = wb[wb.MDepth == depthmin]['dz_ft'].values[0]
        else:
            # Use last col if survey not deep enough
            dx = wb.dx_ft.values[-1]
            dy = wb.dy_ft.values[-1]
            dz = wb.dz_ft.values[-1]

        return (dx, dy, dz)

    def open_lwt(self, file):
        """
        Parses american format lwt (CWR export) sonar file into a DataFrame.

        Parameters
        ----------
        file : str
            Path to the .lwt sonar file.

        Returns
        -------
        df : pandas.DataFrame
            DataFrame parsed from lwt.
        """
        with open(file, 'r') as f:
            text = f.read()

        lines = text.split('\n')

        dpth_i = [i for i, x in enumerate(lines) if 'DEPTH' in x] + [None]

        # Empty pandas.DataFrame
        df = pd.DataFrame()

        for i0, i1 in zip(dpth_i[:-1], dpth_i[1:]):

            tlines = lines[i0:i1]

            # Drop Junk Lines
            tlines = [x.replace('\x0c', '') for x in tlines if 'Page' not in x]
            tlines = [x for x in tlines if len(x) > 0]

            # Parse section header info
            dpth, tilt, rng, vos = [x.strip().split(' ')[0] for x in tlines[0].split(':')][1:]

            # Parse table data into a dataframe for the tilt/depth/range section
            a = []
            for line in tlines[2:]:
                a.append([float(x) for x in line.strip().split(' ') if len(x) > 0])
            tdf = pd.DataFrame(a)

            # Rename columns
            tdf.columns = ['azi', '0.0', '2.8', '5.6', '8.4', '11.3', '14.1', '16.9', '19.7']  # 360 / shot len 2.8125

            # Set header data to temp table
            tdf['depth'] = dpth
            tdf['tilt'] = tilt
            tdf['range'] = rng
            tdf['vos'] = vos

            # Reorder Columns
            tdf = tdf[
                ['depth', 'tilt', 'range', 'vos', 'azi', '0.0', '2.8', '5.6', '8.4', '11.3', '14.1', '16.9', '19.7']]

            # Add to survey wide table
            df = pd.concat([df, tdf], ignore_index=True)

        # self.df = df
        df.rename(columns={'azi': 'cAzi'}, inplace=True)

        for col in df.columns:
            df[col] = df[col].astype(float)

        return df

    def getMagDev(self, long, lat, alt=0, year=2024):
        """
        Computes magnetic declination at a given location.
        Using pyCRGI : https://github.com/pleiszenburg/pyCRGI?tab=readme-ov-file

        D: declination (+ve east) [degree]
        I: inclination (+ve down) [degree]
        H: horizontal intensity [nT]
        X: north component [nT]
        Y: east component [nT]
        Z: vertical component (+ve down) [nT]
        F: total intensity [nT]

        Parameters
        ----------
        long : float
            Longitude in decimal degrees.
        lat : float
            Latitude in decimal degrees.
        alt : float, default=0
            Altitude in km.
        year : int, default=2024
            Year of the survey.

        Returns
        -------
        d : float
            Magnetic declination (+ve east).
        """

        # Get declination and store and return
        d, i, h, x, y, z, f = get_value(lat, long, alt, year)
        self.magdec = d
        return d

    def lwt_df_to_dfT(self, df):
        """
        Converts unformated pandas.DataFrame to one row per sonar shot.

        Parameters
        ----------
        df : pandas.DataFrame
            Sonar data parsed from LWT file.
        Returns
        -------
        dfT : pandas.DataFrame
            Reformatted long-form DataFrame where each row represents a single sonar shot,
            including attributes such as depth, range, azimuth (cAzi), and tilt.
        """

        # Seperate offsets to rows
        dfT = pd.DataFrame()

        for index, row in df.iterrows():
            for i, aziplus in enumerate(['0.0', '2.8', '5.6', '8.4', '11.3', '14.1', '16.9', '19.7']):
                d = {'depth': [row['depth']],
                     'tilt': [row['tilt']],
                     'range': [row['range']],
                     'vos': [row['vos']],
                     'cAzi': [row['cAzi'] + (i * 2.8125)],
                     'r': [row[aziplus]]}
                t = pd.DataFrame(d)
                dfT = pd.concat([dfT, t], ignore_index=True)

        for col in list(dfT):
            dfT[col] = dfT[col].astype(float)

        return dfT
    
    def read_lwt_dat_export(self, path):
        """
        Reads a .dat export from CavView II and reformats into LWT
        
        Parameters
        ----------
        path : str
            File path to the `.dat` file exported from CavView II.

        Returns
        -------
        dfT : pandas.DataFrame
            Reformatted long-form DataFrame where each row represents a single sonar shot,
            including attributes such as depth, range, azimuth (cAzi), and tilt.
        """
        
        # Open File
        with open(path, 'r') as f:
            lines = f.readlines()
        
        # Find Depth lines
        dindex = [i for i,x in enumerate(lines) if 'Depth:' in x]
        dindex = dindex + [None]
        
        # Empty container to add data to
        dfT = pd.DataFrame()
        
        # Loop through lwt sub tables
        for start,end in zip(dindex[:-1], dindex[1:]):
            # Grab and clean subtable
            s = lines[start:end]
            s = [x.replace('\n','') for x in s]
            s = [x for x in s if len(x) > 0]
            
            # Depth
            depth = float(s[0].split(' ')[1])
            
            # Grab and clean subtable
            data = [x.split(' ') for x in s[2:]]
            data = [[y for y in x if y != ''] for x in data]
            data = np.array(data)
            data = data.astype(float)
            
            # Calculate radii steps to handle mutliple outputs
            basedeg = data[:,0]
            data = data[:,1:]
            n = data.shape[1]
            delta = basedeg[1] - basedeg[0]
            step = delta / n

            # Build pandas.DataFrame
            r = data.flatten()
            deg = np.arange(0,360,step)
            n = deg.shape[0]
            
            temp = pd.DataFrame({'depth':[depth]*n,
                                 'tilt':[0]*n,
                                 'range':[np.nan]*n,
                                 'vos':[np.nan]*n,
                                 'cAzi':deg,
                                 'r':r,
                                 })
            
            # Concat to master
            dfT = pd.concat([dfT, temp], ignore_index=True)
        
        # Convert to float
        dfT = dfT.astype(float)

        return dfT
    
    def dfT_to_xyz_delta_points(self, dfT):
        """
        Converts dfT (one row per sonar shot with math inc and azi) and calculates delta XYZ coordinates.

        Parameters
        ----------
        dfT : pandas.DataFrame
            Reformatted long-form DataFrame where each row represents a single sonar shot,
            including attributes such as depth, range, azimuth (cAzi), and tilt.

        Returns
        -------
        xyz : pandas.DataFrame
            Formatted pandas DataFrame with calculated dx, dy, dz (delta position) information.

        Notes
        -----
        Conversion from spherical to Cartesian coordinates using:

            r : radius  
            θ : azimuth  
            φ : inclination  

            x = r * sin(φ) * cos(θ)  
            y = r * sin(φ) * sin(θ)  
            z = r * cos(φ)
        """
        ## Magnetic Declination
        dfT['cAziN'] = dfT['cAzi'] + self.magdec
        dfT['mAzi'] = dfT['cAziN'].apply(self.cAzi2mAzi)
        dfT['mInc'] = 90 - dfT['tilt']

        xyz = dfT.copy()

        ## Add check for 0 and 180 mInc, might be not tilt and mInc
        if xyz[xyz.mInc == 0].shape[0] + xyz[xyz.mInc == 180].shape[0] > 0:
            warnings.warn('Double check mInc as it may be "tilt" because a 0 or 180 shot was detected.')

        # Math done in spherical math azimuth - Use cAzi2mAzi to convert from compass azimuth
        xyz['dx'] = xyz['r'] * np.sin(np.deg2rad(xyz['mInc'])) * np.cos(np.deg2rad(xyz['mAzi']))
        xyz['dy'] = xyz['r'] * np.sin(np.deg2rad(xyz['mInc'])) * np.sin(np.deg2rad(xyz['mAzi']))
        xyz['dz'] = xyz['r'] * np.cos(np.deg2rad(xyz['mInc']))

        return xyz

    def utm_epsg_code_at_point(self, surveyxyz, crs='EPSG:4326'):
        """
        Determines the WGS 84 UTM EPSG code for a given geographic point and 
        transforms the survey coordinates to that UTM zone.

        Parameters
        ----------
        surveyxyz : tuple
            The (x, y, z) of the survey datum in the input CRS.
        crs : str, default='EPSG:4326'
            Coordinate reference system of the input surveyxyz point.

        Notes
        ----
        self.zonecrs : str
            EPSG code string of the identified UTM zone.
        self.pntgdf : GeoDataFrame
            Input point as a GeoDataFrame in its original CRS.
        self.pntgdfUTM : GeoDataFrame
            Input point reprojected to the UTM CRS.
        self.surveyxyzUTM : tuple
            (x, y, z) coordinates transformed into UTM.
        """

        # Wrap the input coordinates as a GeoDataFrame
        pntgdf = gpd.GeoDataFrame({'geometry': [Point(surveyxyz)]},
                                  crs=crs)

        # Make sure CRS is consistant
        if pntgdf.crs != self.utmzones.crs:
            pntgdf.to_crs(self.utmzones.crs, inplace=True)

        self.pntgdf = pntgdf

        # Find appropriate UTM EPSG Code
        zone = self.utmzones[self.utmzones.intersects(pntgdf.unary_union)].ZONE.values[0]
        zone = str(zone)

        while len(zone) < 2:
            zone = '0' + zone

        zonecrs = 'EPSG:326' + zone

        self.zonecrs = zonecrs

        # Convert input point to correct UTM CRS
        self.pntgdfUTM = pntgdf.to_crs(self.zonecrs)

        # Extract transformed UTM coords and convert Z to meters
        self.surveyxyzUTM = list(self.pntgdfUTM.geometry.values[0].coords)[0]
        self.surveyxyzUTM = (self.surveyxyzUTM[0],
                             self.surveyxyzUTM[1],
                             self.surveyxyzUTM[2] / 3.28084)

    def generate_gdf(self, xyz, wb_delta=(0, 0, 0)):
        """
        Converts a DataFrame of delta XYZ values into two GeoDataFrames: one in UTM and one in the specified output CRS.

        Parameters
        ----------
        xyz : pandas.DataFrame
            DataFrame with delta 'dx', 'dy', 'dz', and 'depth' columns for each sonar point.
        wb_delta : tuple of float, optional
            (dx, dy, dz) correction from well deviation survey, default is (0, 0, 0).

        Returns
        -------
        xyz_gdf : geopandas.GeoDataFrame
            Point cloud in the UTM CRS defined by the survey location.
        xyz_gdfCRS : geopandas.GeoDataFrame
            Point cloud in the specified output CRS.
        """
        # Extract survey datum from the GeoDataFrame
        x, y, z = tuple(self.pntgdf.geometry[0].coords)[0]

        # Apply wellbore deviation correction
        wb_dx, wb_dy, wb_dz = wb_delta
        xyz['dx'] += wb_dx
        xyz['dy'] += wb_dy
        xyz['dz'] += wb_dz

        # Convert from feet to meters if working in metric mode
        if not self.metric:
            z /= 3.28084
            for col in ['dx', 'dy', 'dz']:
                xyz[col] = xyz[col].astype(float) / 3.28084

        # Compute full UTM coordinates
        xyz['x'] = self.surveyxyzUTM[0] + xyz['dx']
        xyz['y'] = self.surveyxyzUTM[1] + xyz['dy']
        xyz['z'] = self.surveyxyzUTM[2] - ((xyz['depth'] / 3.28084) - xyz['dz'])

        # Create 3D shapely Points (loop version, explicit)
        xyz['geometry'] = ''
        for index, row in xyz.iterrows():
            xyz.at[index, 'geometry'] = Point(row['x'], row['y'], row['z'])

        # Create GeoDataFrame in UTM
        xyz_gdf = gpd.GeoDataFrame(xyz, crs=self.zonecrs)

        # Reproject to target CRS
        xyz_gdfCRS = xyz_gdf.copy().to_crs(self.crs)

        # Fix z-axis units if metric is False (since CRS transform doesn't touch z)
        if not self.metric:
            for index, row in xyz_gdfCRS.iterrows():
                x, y, z = list(row['geometry'].coords)[0]
                z *= 3.28084
                xyz_gdfCRS.at[index, 'geometry'] = Point(x, y, z)
                xyz_gdfCRS.at[index, 'x'] = x
                xyz_gdfCRS.at[index, 'y'] = y
                xyz_gdfCRS.at[index, 'z'] = z

        # Store attributes for reuse
        self.xyz_gdf = xyz_gdf
        self.xyz_gdfCRS = xyz_gdfCRS

        return xyz_gdf, xyz_gdfCRS


    def cAzi2mAzi(self, a):
        """
        Convert compass azimuth to math spherical azimuth.

        Parameters
        ----------
        a : float
            Compass azimuth in degrees

        Returns
        -------
        d : float
            Math spherical azimuth
        """

        # Return NaN if input is not a number
        if pd.isna(a):
            return np.nan
        
        # Convert compass to math azimuth: math = 90 - compass
        d = (90 - a) % 360

        return d

    def degmmss_to_dec_deg(x):
        """
        Converts DMS (degrees minutes seconds) strings to decimal degrees.
        
        Parameters
        ----------
        x : string
            Degrees minutes seconds in format : 93°24'50.61"W

        Returns
        -------
        dd : float
            Decimal degrees
        """
        # confirm " isn't ''
        x = x.replace("''", '"')

        # extract DMS
        d = float(x.split('°')[0])
        m = float(x.split('°')[1].split("'")[0])
        s = float(x.split("'")[1].split('"')[0])

        # Convert do decimal degrees
        dd = d + (m / 60) + ((s / 60) / 60)

        if 'W' in x:
            dd = dd * -1
        elif "S" in x:
            dd = dd * -1

        return dd

class SonarPlots:
    """
    Python Class SonarPlots

    A collection of methods for processing, visualizing, and analyzing sonar cavern surveys.

    ## Example ------------------------------------------------------------------------------------

    """
    
    def ckdnearest(gdfA, gdfB, gdfB_cols=['Place']):
        """
        Find the nearest neighbors between two point cloud GeoDataFrames
        and return a GeoDataFrame of the line between nearest points.
        
        For estimating the web thickness.

        Parameters
        ----------
        gdfA: geopandas.GeoDataFrame
            Source point cloud.
        gdfB: geopandas.GeoDataFrame
            Target point cloud.
        gdfB_cols: list, optional
            Columns to return from gdfB. Default is ['Place'].

        Returns
        -------
        line: geopandas.GeoDataFrame
            Line between nearest points in gdfA and gdfB
        """
        # Remove known error sources and reset indexing
        gdfA = gdfA.copy()
        gdfB = gdfB.copy()
        gdfA.reset_index(drop=True, inplace=True)
        gdfB.reset_index(drop=True, inplace=True)
        
        # Drop leftover columns if re-run
        if 'geometry1' in gdfA.columns:
            gdfA.drop(columns=['geometry1'], inplace=True)
        if 'geometry1' in gdfB.columns:
            gdfB.drop(columns=['geometry1'], inplace=True)

        # Add temporary second geometry column
        gdfB['geometry1'] = gdfB.geometry
        gdfB_cols.append('geometry1')

        # Extract point coordinates
        A = np.concatenate([np.array(geom.coords) for geom in gdfA.geometry.to_list()])
        B = [np.array(geom.coords) for geom in gdfB.geometry.to_list()]

        # Build mapping from each point to its row index
        B_ix = tuple(itertools.chain.from_iterable(
            [itertools.repeat(i, x) for i, x in enumerate(list(map(len, B)))]))
        B = np.concatenate(B)

        # Build cKDTree for fast nearest neighbor search
        ckd_tree = cKDTree(B)
        dist, idx = ckd_tree.query(A, k=1)
        idx = itemgetter(*idx)(B_ix)

        # Combine results: source point, nearest neighbor, and distance
        gdf = pd.concat([
            gdfA,
            gdfB.loc[idx, gdfB_cols].reset_index(drop=True),
            pd.Series(dist, name='dist')
        ], axis=1)

        # Identify the closest match
        i = gdf[gdf.dist == gdf.dist.min()].index.values[0]
        p0, p1 = gdf.at[i, 'geometry'], gdf.at[i, 'geometry1']

        # Handle MultiPoint geometries if needed
        if isinstance(p1, Point):
            pass
        elif len(list(p1)) > 1:
            p1 = p1[0]

        # Create line between matched points
        line = gpd.GeoDataFrame({
            'geometry': [LineString([p0, p1])],
            'distance': [gdf.dist.min()]
        }, crs=gdfA.crs)

        return line


    def addXY2RadiusAziM(gdf):
        """
        Adds radius and azimuth (math) to a pandas.DataFrame or geopandas.GeoDataFrame.

        Parameters
        ----------
        gdf : pandas.DataFrame or geopandas.GeoDataFrame
            Must contain 'X' and 'Y', or 'dx' and 'dy'.

        Returns
        -------
        gdf : pandas.DataFrame
            The input DataFrame with added columns:
                - 'mAzi': mathematical azimuth
                - 'cAzi': compass azimuth
                - 'radiusXY': 2D radius from origin in XY plane
                - 'radiusXYZ': 3D radius from origin
        """
        
        if 'X' not in gdf.columns:
            rename = {'dx': 'X',
                      'dy': 'Y',
                      'dz': 'Z'}
            gdf.rename(columns=rename, inplace=True)
        # Check inputs
        if 'X' not in gdf.columns:
            raise Exception('X column not in input')
        if 'Y' not in gdf.columns:
            raise Exception('Y column not in input')

        # math azimuth
        gdf['mAzi'] = np.degrees(np.arctan2(gdf['Y'], gdf['X']))
        gdf['mAzi'] = gdf['mAzi'].where(gdf['mAzi'] >= 0, other=gdf['mAzi'] + 360)

        # Compass Azimuth
        gdf['cAzi'] = gdf['mAzi'].apply(SonarPlots.cAzi_2from_mAzi)

        # Radius
        gdf['radiusXY'] = np.sqrt(gdf['X'] ** 2 + gdf['Y'] ** 2)
        gdf['radiusXYZ'] = np.sqrt(gdf['X'] ** 2 + gdf['Y'] ** 2 + gdf['Z'] ** 2)

        return gdf

    def cAzi_2from_mAzi(a):
        """
        Converts between compass azimuth and mathematical azimuth.

        Parameters
        ----------
        a : float
            Math sperical degrees or compass azimuth degrees

        Returns
        -------
        d : float
            Converted azimuth in degrees.
        """
        if pd.isna(a):
            return np.nan

        # Perform conversion
        d = (90 - a) % 360

        return d

    def xyzexterior(gdf, overridebuff=None):
        """
        Creates an exterior shapely.geometry.Polygon from and geopandas.GeoDataFrame xyz point cloud.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            xyz point cloud.
        overridebuff : float, optional
            Adds a buffer rather than optimal distance. Use if shape is irregular.
        
        Returns
        -------
        buff : shapely.geometry.Polygon
            Exterior polygon of input points.
        """
        # Add radius info
        gdf = SonarPlots.addXY2RadiusAziM(gdf)

        # Unit conversion
        # for col in ['radiusXY','radiusXYZ']:
        #    gdf[col] = gdf[col] / unitconversion

        # Find min buffer
        rmax = gdf['radiusXY'].max()
        buffMin = rmax * np.tan(np.deg2rad(2.8125))

        # Override buffer?
        if overridebuff:
            print('switch to', overridebuff)
            buffMin = overridebuff

        # Build Buffer Geometry
        buff = gdf.buffer(buffMin * 2, join_style=3)
        buff = buff.unary_union.buffer(buffMin * -2, join_style=3)
        buff = SonarPlots.remove_interiors(buff)

        return buff


    def remove_interiors(poly):
        """
        Remove interior holes from a polygon

        Parameters
        ----------
        poly : shapely.geometry.Polygon
            Input shapely Polygon.

        Returns
        ---------
        poly : shapely.geometry.Polygon
            Polygon without any interior holes.
        """
        if isinstance(poly, MultiPolygon):
            return poly.geoms[0]

        elif poly.interiors:
            return Polygon(list(poly.exterior.coords))

        else:
            return poly

    def volume_df_from_stl(polydatamesh, dstep=10, correction=True):
        """
        Computes the volume profile of a 3D cavern mesh along the Z-axis.

        Parameters
        ----------
        polydatamesh : pyvista.core.pointset.PolyData
            Mesh of cavern.
        dstep : int, optional
            Depth step interval for slicing along the Z-axis. Default is 10.
        correct : bool, optional
            If True, scales computed volumes so the final value matches the mesh's total volume.

        Returns
        -------
        df : pandas.DataFrame
            df of depth and cumulative volume in cubic units.

        RETURNS
        df: pandas.DataFrame, Dataframe of depth and cumulative volume in cubic units of mesh.

        v2: updated threshold and grid creation
        """
        xmin, xmax, ymin, ymax, zmin, zmax = polydatamesh.bounds

        # Make grid
        xx = np.arange(int(xmin - 0.5), int(xmax + 2))
        yy = np.arange(int(ymin - 0.5), int(ymax + 2))
        zz = np.arange(int(zmin - 0.5), int(zmax + 2))
        dataset = pv.RectilinearGrid(xx, yy, zz)
        dataset.compute_implicit_distance(polydatamesh, inplace=True)

        # Select mesh interior
        inner = dataset.threshold(0.51, scalars="implicit_distance", invert=True)

        # Compute volumes along Z
        depths = np.arange(zmin - dstep, zmax + dstep, dstep)
        volumes = []
        for d in depths:
            clip_bounds = [xmin - 5, xmax + 5, ymin - 5, ymax + 5, zmin - 5, d]
            clipped = inner.clip_box(clip_bounds)
            volumes.append(clipped.volume)
            print('.', end='')

        df = pd.DataFrame({'depth': depths, 'volume': volumes})

        # Optional correction to match total mesh volume
        if correction:
            correction_factor = polydatamesh.volume / df['volume'].max()
            df['volume'] = (df['volume'] * correction_factor).astype(int)

        return df

    def old_volume_df_from_stl(polydatamesh, dstep=10):
        """
        Provides a z-axis volume profile of a cavern in the units**3 of the cavern mesh (usually m of ft).

        Parameters
        ----------
        polydatamesh : pyvista.core.pointset.PolyData
            3D mesh of cavern geometry.
        dstep : int, optional
            Vertical interval (step size) to slice and compute volume, default is 10.

        Returns
        -------
        df : pandas.DataFrame
            Dataframe of depth and cumulative volume in cubic units of mesh.
        """
        # convert mesh to grid
        grid = polydatamesh.cast_to_unstructured_grid()

        # Get mesh bounds
        xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
        dx = int(xmax - xmin + 1)
        dy = int(ymax - ymin + 1)
        dz = int(zmax - zmin + 1)

        # Make voxel 3D grid
        grid = pv.ImageData()
        grid.dimensions = np.array([dx, dy, dz])
        grid.spacing = (1, 1, 1)
        grid.origin = (xmin, ymin, zmin)

        grid.compute_implicit_distance(polydatamesh, inplace=True)
        grid = grid.threshold(0.0, scalars="implicit_distance", invert=True)

        # Compute Volumes
        depths = np.arange(zmin, zmax, dstep)
        volumes = []
        for d in depths:
            clip_bounds = [xmin, xmax, ymin, ymax, zmin, d]
            clipped = grid.clip_box(clip_bounds)
            vol = clipped.volume
            volumes.append(vol)
            print('.', end='')

        # df of volume at each depth
        df = pd.DataFrame({'depth': depths, 'volume': volumes})

        return df

    def half_cav_rgba(path, cav=None, xy0=None):
        """
        Generates a clipped image of a cavern mesh from an STL file.

        Loads a 3D mesh of a cavern, slices it vertically at a given (x0, y0) point (or automatically at the highest point),
        and returns a cropped RGBA image.

        Parameters
        ----------
        path : str
            File path to STL mesh.
        cav : pyvista.PolyData, optional
            Already-loaded PyVista mesh of the cavern. If not provided, it will be read from `path`.
        xy0 : tuple of float, optional
            Tuple (x0, y0) defining the vertical clip location. If not provided, the mesh's topmost point is used.

        Returns
        -------
        crop_img : np.array
            RGBA NumPy image array of the clipped mesh, cropped to remove white margins.
        """
        # Open if needed
        if cav is None:
            cav = pv.read(path)
        if xy0:
            x0, y0 = xy0
        else:
            # Find highest Z point in mesh
            m = mesh.Mesh.from_file(path)
            z1 = m.points[:, 2::3].min()
            z0 = m.points[:, 2::3].max()
            top = m.points[m.points[:, 2] == z0]
            x0 = top[:, 0::3].mean()
            y0 = top[:, 1::3].mean()

        xmin, xmax, ymin, ymax, zmin, zmax = cav.bounds
        clip_bounds = [xmin, x0, ymin, ymax, zmin, zmax]
        clipped = cav.clip_box(clip_bounds)
        pl = pv.Plotter(off_screen=True)
        pl.add_mesh(clipped, color='#328da8')
        pl.camera_position = 'xz'
        img = pl.screenshot()  # doesn't open display like other method

        # Compute alpha channel from RGB values
        img_alpha = np.mean(img[:, :, :3], axis=2)  # Compute average intensity as alpha

        # Threshold alpha channel to create binary mask (optional)
        ma = img_alpha == 255
        m = img_alpha < 255
        img_alpha[ma] = 0
        img_alpha[m] = 255

        # Expand dimensions of img_alpha to match img for concatenation
        img_alpha = np.expand_dims(img_alpha, axis=2)

        # Concatenate img and img_alpha along the third axis to create RGBA image
        rgba_img = np.concatenate((img, img_alpha), axis=2).astype(int)
        pl.close()

        ## Crop whitespace from left and right sides
        s0 = rgba_img.mean(axis=0).sum(axis=1) != 765

        crop_img = rgba_img[:, s0]

        return crop_img

    def read_inv_file(path):
        """
        Reads a Sonarwire inventory file and parses it into a pandas DataFrame.

        Parameters
        ----------
        path : str
            File path to inventory file.

        Returns
        -------
        df : pandas.DataFrame
            DataFrame of inventory depths and volumes in columns:
            ['Depth', 'BBL_ft', 'BBL']
        """
        with open(path, 'r') as f:
            lines = f.readlines()

        data = [[y.replace('\n', '') for y in x.split(' ') if len(y) > 0] for x in lines]
        data = [x for x in data if len(x) == 3]
        data = [x for x in data if x[0].isnumeric()]
        df = pd.DataFrame(data)
        df.columns = ['Depth', 'BBL_ft', 'BBL']
        df = df.astype(float)

        return df


    def parse_depth_chart(df, dstep=10):
        """
        Parses a depth-volume DataFrame on regular dstep intervals
        for cumulative volume plotting.

        * Missing z component; change to z space outside function to be the same as SonarPlots.volume_df_from_stl output.

        Parameters
        ----------
        df : pandas.DataFrame
            Ouput from read_inv_file()
        dstep : int, optional)
            Step interval for sampling, default 10

        Returns
        -------
        t : pandas.DataFrame
            DataFrame sampled on dstep spacing.
        """

        dmin = df.Depth.min()
        dmax = df.Depth.max()
        a = np.arange(dmin, dmax + dstep, dstep)
        t = df[df.Depth.isin(a)][['Depth', 'BBL']].copy()

        return t


class SonarDXF:
    """
    Class for processing and interpolating sonar 3D line data.
    """
    def alpha_infill_points(geom, alpha=30):
        """
        Generates infill points along LineString based on spacing threshold (alpha).
        These are to attempt to fill voids from meshing.

        Parameters
        ----------
        geom : shapely.geometry.Linestring
            3D linestring object.
        alpha : float, optional
            Spacing threshold between infill points, default is 30.

        Returns
        -------
        points : list of shapely.geometry.Points
            List of generated infill points along the LineString.
        """
        xyz = list(geom.coords)
        points = []
        
        for xyz0, xyz1 in zip(xyz[:-1], xyz[1:]):

            ## Infill point count
            x0, y0, z0 = xyz0
            x1, y1, z1 = xyz1
            d = ((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2) ** 0.5
            n = d // alpha

            if n == 0:
                continue

            # Make infill points
            xi = np.linspace(xyz0, xyz1, int(n + 2))
            t = xi[1:-1]
            t = [Point(p) for p in t]
            points = points + t

        return points

    def infill_points(geom, n=3):
        """
        Generates a fixed numebrer (n) of infill points along LineString.
        These are to attempt to fill voids from meshing.

        PArameters
        ----------
        geom : shapely.geometry.Linestring
            3D linestring object.
        alpha : int, optional
            Number of points to insert between each segment, default 3.

        Returns
        -------
        points : list of shapely.geometry.Points
            List of generated infill points along the LineString.
        """
        xyz = list(geom.coords)
        points = []

        for xyz0, xyz1 in zip(xyz[:-1], xyz[1:]):
            xi = np.linspace(xyz0, xyz1, n + 2)
            t = xi[1:-1]
            t = [Point(p) for p in t]
            points += t

        return points

    def make_infill_df(lines, alpha=20):
        """
        Generates a DataFrame of infill points from line geometries

        Parameters
        ----------
        lines : pd.DataFrame or geopandas.GeoDataFrame
            DataFrame with 'geometry' column containing 3D LineString objects.
        alpha : float, optional
            Spacing threshold between infill points, default is 320.
        
        Returns
        -------
        infill : pandas.DataFrame
            DataFrame with infill point geometries.
        """
        points = []

        # loop through each line and generate infill points
        for index, row in lines.iterrows():
            t = SonarDXF.alpha_infill_points(row['geometry'], alpha=alpha)
            points = points + t

        # convert to DataFrame
        infill = pd.DataFrame()
        infill['geometry'] = points

        # Extract xyz from geometry
        for index, row in infill.iterrows():
            x, y, z = list(row['geometry'].coords)[0]
            infill.at[index, 'x'] = x
            infill.at[index, 'y'] = y
            infill.at[index, 'z'] = z

        return infill

    def process_lines(xyz_gdf):
        """
        Converts a GeoDataFrame of 3D points into horizontal and vertical LineString segments
        
        Parameters
        ----------
        xyz_gdf : pandas.DataFrame or geopandas.GeoDataFrame
            DataFrame of points with columns: 'x', 'y', 'z', 'depth', 'tilt', and 'cAziN'.
        
        Returns
        -------
        lines : pandas.DataFrame
            DataFrame with a 'geometry' column containing 3D LineString objects.
        """
        lines = pd.DataFrame()

        # Horizontal lines
        for (depth, inc), group in xyz_gdf.groupby(['depth', 'tilt']):
            group.sort_values('cAziN', inplace=True)
            pnts = group[['x', 'y', 'z']].values
            pnts = np.vstack([pnts, pnts[0]])
            geom = LineString(pnts)
            t = pd.DataFrame({'depth': [depth], 'inc': [inc],})
            t['geometry'] = geom
            t['z_m'] = group.z.mean()
            lines = pd.concat([lines, t], ignore_index=True)

        # Verticle lines : tilt shots
        depths = xyz_gdf[xyz_gdf.tilt != 0].depth.unique()
        for depth in depths:
            group = xyz_gdf[xyz_gdf.depth == depth].copy()
            group['_azi'] = group['cAziN'].where(group['cAziN'] < 180, other=group['cAziN'] - 180)

            for azi, g2 in group.groupby('_azi'):
                ## Sort points by inc, one side then the other
                g20 = g2[g2.cAziN < 180].sort_values('tilt')
                g21 = g2[(g2.cAziN >= 180) & (g2.index.isin(g20.index.tolist()))].sort_values('tilt', ascending=False)

                g2 = pd.concat([g2[g2.cAziN < 180].sort_values('tilt'),
                                g2[g2.cAziN >= 180].sort_values('tilt', ascending=False)],
                               ignore_index=True)

                pnts = g2[['x', 'y', 'z']].values
                t = pd.DataFrame({'depth': [depth],'inc': [np.nan]})
                t['cAziN'] = azi
                t['geometry'] = LineString(pnts)
                t['z_m'] = group.z.mean()

                lines = pd.concat([lines, t], ignore_index=True)

        # Verticle lines : flat shots (inc == 0)
        for azi, group in xyz_gdf[xyz_gdf.tilt == 0].groupby('cAziN'):
            group.sort_values('z', ascending=False, inplace=True)
            pnts = group[['x', 'y', 'z']].values

            t = pd.DataFrame({'depth': [depth],'inc': [np.nan]})
            t['geometry'] = LineString(pnts)
            t['cAziN'] = azi
            t['z_m'] = group.z.mean()

            lines = pd.concat([lines, t], ignore_index=True)

        return lines

    def extract_points(lines):
        """
        Extracts 3D coordinate points from LineString geometries in a GeoDataFrame.

        Parameters
        ----------
        lines : pandas.DataFrame or geopandas.GeoDataFrame
            DataFrame with 'geometry' column of LineString objects

        Returns
        -------
        linepoints : list of tuple
            Flattened list of (x, y, z) coordinate tuples.
        """

        nested_list = lines.geometry.apply(lambda x: list(x.coords)).tolist()
        linepoints = [item for sublist in nested_list for item in sublist]

        return linepoints

class SonarPyVista:
    """
    Class for creating PyVista-based visualizations from sonar-derived cavern data.
    """
    def __init__(self):
        self.__version__ = '0.0.0a'

    def wireframe_from_cavlines2(lines):
        """
        Builds a PyVista wireframe mesh from cavern survey line geometries.

        Parameters
        ----------
        lines : pandas.DataFrame or geopandas.GeoDataFrame
            DataFrame with 'geometry' column with 3D LineString objects.

        Returns
        -------
        mesh : pyvista.PolyData
            Merged wireframe mesh.
        """

        data = []
        for index, row in lines.iterrows():
            poly_line = pv.MultipleLines(points=list(row['geometry'].coords))
            data.append(poly_line)

        mesh = pv.merge(data)

        return mesh


def cwt2lwt(path, exepath):
    """
    Automates the Sonarwire Postprocessing software to convert a .CWR file into a longwall table (.LWT) report.

    Parameters
    ----------
    path : str
        Filepath to CWR file.
    exepath: str
        Filepath to sonarpostprocessing.exe.
    """
    default_cwd = os.getcwd() + '/default.sdb'
    if os.path.exists(default_cwd):
        os.remove(default_cwd)

    ## Open app
    process = subprocess.Popen([exepath])
    time.sleep(1.25)

    ## Load File
    hwnd = win32gui.FindWindow(None, 'Sonarwire Postprocessing')
    pyautogui.press("alt", interval=0.1)
    win32gui.SetForegroundWindow(hwnd)
    time.sleep(2)
    x0, y0, x1, y1 = win32gui.GetWindowRect(hwnd)
    vieweditfield = pyautogui.locateOnScreen('cav.png')
    vef_center = pyautogui.center(vieweditfield)
    pyautogui.moveTo(vef_center)
    pyautogui.click()

    ## Choose location
    time.sleep(1.5)
    pyautogui.keyDown('ctrl')
    pyautogui.press("l", interval=0.1)
    pyautogui.keyUp('ctrl')

    pyautogui.write('/'.join(path.split('/')[:-1]) + '/')
    pyautogui.press("enter", interval=0.1)
    pyautogui.press("tab", interval=0.1)
    pyautogui.press("tab", interval=0.1)
    pyautogui.press("tab", interval=0.1)
    pyautogui.press("tab", interval=0.1)
    pyautogui.press("tab", interval=0.1)
    pyautogui.press("tab", interval=0.1)
    time.sleep(0.5)
    pyautogui.write(path.split('/')[-1])

    time.sleep(0.2)
    pyautogui.press("enter", interval=0.1)
    time.sleep(1.0)
    if pyautogui.getActiveWindow().title != 'Survey File Parameters':
        pyautogui.press("enter", interval=0.1)  # seems to have fixed that hang up
        time.sleep(2)

    ## Accept default survey file parameters
    time.sleep(0.5)
    pyautogui.press("tab", interval=0.2)
    pyautogui.press("tab", interval=0.2)
    pyautogui.press("enter", interval=0.2)
    time.sleep(1.5)

    ## Check for error
    if win32gui.FindWindow(None, 'Error') != 0:
        print(win32gui.FindWindow(None, 'Error'), win32gui.FindWindow(None, 'Error') == 0)
        return None

    # Close Sonar view window
    time.sleep(0.5)
    pyautogui.press("alt", interval=0.1)
    pyautogui.press("f", interval=0.1)
    pyautogui.press("c", interval=0.1)
    time.sleep(0.5)

    # Set up report export
    pyautogui.press("alt", interval=0.1)
    pyautogui.press("u", interval=0.1)
    pyautogui.press("r", interval=0.1)
    time.sleep(1)
    pyautogui.press("down", interval=0.1)

    x0, y0, x1, y1 = win32gui.GetWindowRect(win32gui.GetForegroundWindow())
    reporttype = pyautogui.locateOnScreen('swt.png')
    rt_center = pyautogui.center(reporttype)
    pyautogui.moveTo(rt_center)
    pyautogui.click()

    time.sleep(0.2)
    pyautogui.press("down", interval=0.2)
    pyautogui.press("down", interval=0.2)
    pyautogui.press("down", interval=0.2)
    pyautogui.press("down", interval=0.2)
    pyautogui.press("down", interval=0.2)
    pyautogui.press("enter", interval=0.2)

    pyautogui.press("tab", interval=0.2)
    pyautogui.press("enter", interval=0.2)
    time.sleep(1.5)

    # Sometimes the report window opens — this is still a problem
    reportwindow = [x.title for x in pyautogui.getAllWindows() if 'Sonarwire Report Viewer:' in x.title]
    if len(reportwindow) > 0:
        time.sleep(1)
        print('^')

        hwnd = win32gui.FindWindow(None, reportwindow[0])
        pyautogui.press("alt", interval=0.1)
        win32gui.SetForegroundWindow(hwnd)
        time.sleep(0.25)

        pyautogui.keyDown('alt')
        pyautogui.press("F4", interval=0.1)
        pyautogui.keyUp('alt')
        time.sleep(1)

    # Exit Report
    pyautogui.press("tab", interval=0.1)
    pyautogui.press("enter", interval=0.1)

    # Exit Program
    pyautogui.press("alt", interval=0.1)
    pyautogui.press("f", interval=0.1)
    pyautogui.press("c", interval=0.1)
    pyautogui.press("c", interval=0.1)

    time.sleep(0.1)
    pyautogui.press("enter", interval=0.1)
    time.sleep(0.1)
    pyautogui.press("tab", interval=0.1)
    pyautogui.press("enter", interval=0.1)
    time.sleep(0.1)
    pyautogui.press("enter", interval=0.1)
    time.sleep(0.1)
    pyautogui.press("enter", interval=0.1)


def parse_date(date_string):
    """
    Format date to yyyymmdd format.

    Parameters
    ----------
    date_string : str
        String containing date in format 'Jul 7, 2024'
    
    Returns
    -------
    formatted_date : str or None
        Date reformatted to yyyymmdd, None if no valid date is found.
    """
    months = {
        'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
        'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
        'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'
    }
    # Define the regex pattern to match the date string
    pattern = r'\w{3}[ ]{1,2}\d{1,2}, \d{4}'

    # Match the date string with the regex pattern
    match = re.findall(pattern, date_string)

    if match:
        # Extract components of the date
        day = match[0].split(' ')[1].replace(',', '')
        month = months[match[0].split(' ')[0]]
        year = match[0][-4:]

        # Format the date as yyyymmdd
        formatted_date = f'{year}{month}{day}'

        return formatted_date
    else:
        return None


def add_rind(array_3d, rind_width):
    """
    Add a rind (border) of zeros to a 3D array.

    Parameters
    ----------
    array_3d : np.ndarray
        3D NumPy array to be padded
    rind_widith : int
        Number of zeros to pad on each side of every axis
    
    Returns
    -------
    padded_array : np.ndarray
        Padded 3D array
    """
    padded_array = np.pad(array_3d, pad_width=rind_width, mode='constant', constant_values=0)

    return padded_array

def remove_rind(array_3d, rind_width):
    """
    Remove rind (border) of zeros from a 3D array
    
    Parameters
    ----------
    array_3d : np.ndarray
        3D NumPy array to remove padding
    rind_width : int
        Number of layers to remove from each side of the array.
    
    Returns
    -------
    array_cleaned : np.ndarray
        Array with padding removed.
    """
    # Calculate the slices to remove the rind
    slices = tuple(slice(rind_width, -rind_width) if dim != 0 else slice(None) for dim in range(array_3d.ndim))

    # Slice the array to remove the rind
    # array_cleaned = array_3d[slices]
    array_cleaned = array_3d[rind_width:-1 * rind_width,
                    rind_width:-1 * rind_width,
                    rind_width:-1 * rind_width]
    return array_cleaned

def mesh_xyz(xyz, iterations=3, n_iter=200, rf=0.1):
    """
    Takes a geopandas.GeoDataFrame of processed xyz date (UTM) and creates
    a surface of the exterior of the cavern.  

    Parameters
    ----------
    xyz : geopandas.GeoDataFrame
        GeoDataFrame containing processed sonar data with columns:
        ['x', 'y', 'z', 'dx', 'dy', 'dz'].
    iterations : int, optional
        Number of dilation/erosion passes, default is 3.
    n_iter : int, optional
        Number of smoothing iterations, default is 200.
    rf : float, optional
        Relaxation factor for smoothing (0 disables smoothing).

    Returns
    -------
    envelope : pyvista.PolyData
        The final triangulated mesh surface of the cavern.
    """
    
    ## Make Array -------------------------------------------------------
    # Try using the val - minimum as an index then refrence later on
    minx, miny, maxx, maxy = xyz.total_bounds
    minz, maxz = xyz.z.min(), xyz.z.max()
    dx = int(maxx - minx + 1)
    dy = int(maxy - miny + 1)
    dz = int(maxz - minz + 1)
    new_variable = np.zeros((dx, dy, dz))
    #dx, dy, dz
    
    # Shot origins
    xyz['x0'] = xyz.x - xyz.dx
    xyz['y0'] = xyz.y - xyz.dy
    xyz['z0'] = xyz.z - xyz.dz
    
    ## Check for shots originating outside the pointcloud (upshots)
    n = sum(xyz['z0'] > maxz + 2) + sum(xyz['z0'] < minz - 2)
    if n > 0:
        raise Exception('n:',n,'minz:', minz, 'maxz:', maxz, 'z0min:',xyz.z0.min(), 'z0max:',xyz.z0.max(),
                        'Shot origin outside of point cloud. Upshot/Downshot survey xyz geodataframes need to be concated with the xyz geodataframes(s) of the rest of the survey for this method')


    for index, row in xyz.iterrows():
        x0 = int(xyz['x0'].values[index] - minx)
        y0 = int(xyz['y0'].values[index] - miny)
        z0 = int(xyz['z0'].values[index] - maxz) * -1

        x1 = int(xyz['x'].values[index] - minx)
        y1 = int(xyz['y'].values[index] - miny)
        z1 = int(xyz['z'].values[index] - maxz) * -1

        # Use Bresenham's Line Algorithm to get the indices of the points along the line
        indices = line_nd((x0, y0, z0),
                          (x1, y1, z1))

        # Set the corresponding points in the array to 1
        new_variable[indices] = 1

    new_variable = add_rind(new_variable, iterations + 1) #testing rind

    struct = scipy.ndimage.generate_binary_structure(3,3)

    dilated = scipy.ndimage.binary_dilation(new_variable, 
                                               iterations=iterations, 
                                               structure=struct)

    eroded = scipy.ndimage.binary_erosion(dilated,
                                               iterations=iterations, 
                                               structure=struct)

    eroded = remove_rind(eroded, iterations + 1) #testing rind

    ## Make Voxel --------------------------------------------------------
    grid = pv.ImageData()
    grid.dimensions = np.array([dx, dy, dz]) + 1 # new_variable.shape
    grid.spacing = (1, 1, 1)                     # Adjust as needed
    grid.origin = (minx, miny, minz)  # bottom SW corner, should be UTM 
    #                                 # but StatePlane would work if you 
    #                                 # adjust the iterations
    grid.cell_data["values"] = eroded[:, :, ::-1].flatten(order="F")
    
    thresh = grid.threshold(0.5)
    surface = thresh.extract_geometry()
    
    # Apply Laplacian smoothing -------------------------------------------
    if rf != 0:
        smoothed_surface = surface.smooth(n_iter=n_iter, relaxation_factor=rf)
    else:
        print('raw')
        smoothed_surface = surface
    
    
    # Fix edges -----------------------------------------------------------
    meshfix = mf.MeshFix(smoothed_surface)

    # Repair also fills holes
    meshfix.repair(verbose=True)

    envelope = meshfix.mesh.clean().triangulate()
    
    return envelope

class ubro:
    """
    Class for exporting and transforming sonar point cloud data into
    Underground Borehole Radar Output (UBRO) format.
    """
    def __init__(self):
        self.__version__ = '0.0.0d'

    def pntUTM_2_ubro_shape(path, xyz_gdfUTM=np.array([]), r=999, header=False):
        """
        Exports a xyz point cloud UTM into a tab delimited format for UBRO.
        
        Parameters
        ----------
        path : str
            File path to the xyz_gdfUTM or a point cloud shapefile in UTM.
        xyz_gdfUTM : geopandas.GeoDataFrame, optional
            UTM-based point cloud. If empty, the function will load from the given file path.
        r : float, optional
            Radius (m) used to generate a 16-spoke radial "star" around the cloud center, default: 999.
        header : bool, optional
            If True, includes column headers in the output CSV (default: False).
        """
        ## Load or use fed file ----------------------------------------------
        if xyz_gdfUTM.shape[0] == 0:
            pnts = gpd.read_file(path)
        else:
            pnts = xyz_gdfUTM.copy()

        if pnts.crs.to_dict()['proj'] != 'utm':
            raise Exception('Not a UTM File', pnts.crs)

        ## Make Origin --------------------------------------------------------
        x0 = (pnts['x'] - (pnts['dx'] / 3.28084)).mean()
        y0 = (pnts['y'] - (pnts['dy'] / 3.28084)).mean()
        z0 = ((pnts['depth'] / 3.28084) + pnts['z']).mean()

        ## Make Star ----------------------------------------------------------
        lines = []
        cazis = [x for x in np.arange(0, 360, 360 / 16)]
        mazis = [SonarPy.cAzi2mAzi(None, a) for a in cazis]

        for angle in mazis:
            x1 = x0 + r * np.cos(np.deg2rad(angle))
            y1 = y0 + r * np.sin(np.deg2rad(angle))
            lines.append(LineString([(x0, y0), (x1, y1)]))

        star = gpd.GeoDataFrame({'cAzi': cazis, 'mAzis': mazis, 'geometry': lines},
                                crs=pnts.crs)
        star.sort_values('cAzi', inplace=True)

        ## Iterrate Through Depths---------------------------------------------
        zmin = pnts.z.min()
        zmax = pnts.z.max()
        depths = range(int(zmin - 1), int(zmax + 1))

        ubroout = pd.DataFrame()
        for z in depths:
            pselect = pnts[(pnts.z > z - 1) & (pnts.z < z + 1)].copy()
            if pselect.shape[0] == 0:
                continue
            r = pselect.r.max() / 3.28084

            exterior = SonarPlots.xyzexterior(pselect,
                                              overridebuff=25)

            c = star.clip(exterior).copy()
            c['range'] = c.geometry.length

            # Missing values (yup that was it)
            for a in [a for a in cazis if a not in c.cAzi.values]:
                c = pd.concat([c, pd.DataFrame({'cAzi': [a],
                                                'range': [0]})],
                              ignore_index=True)

            c['depth'] = z
            ubroout = pd.concat([ubroout, c[['depth', 'cAzi', 'range']]],
                                ignore_index=True)

        ## Transform To UBRO Format -------------------------------------------
        colintmap = {a: i for i, a in enumerate(cazis)}
        ubroout['i'] = ubroout['cAzi'].map(colintmap)
        ubrooutT = pd.DataFrame()
        for i, group in ubroout.groupby('i'):
            # print(i,group.shape[0])
            group['H'] = float(z0) - group['depth']
            group.sort_values('H', inplace=True)
            ubrooutT[f'H{i + 1} [m]'] = group['H'].values
            ubrooutT[f'R{i + 1} [m]'] = group['range'].values

        ## Save out -----------------------------------------------------------
        ubrooutT = ubrooutT.round(3)
        path = path.replace('\\', '/')
        # outpath = path[:path[::-1].index('/')]
        outpath = path.replace('.shp', '_ubroShape.csv')
        ubrooutT.to_csv(outpath, index=False, header=header)

    def ubro_shape2dfT(path, outunit='ft'):
        """
        Reads an UMBRO shape.txt file, infills the azimuth for the 16 sectors.

        Parameters
        ----------
        path : str
            path of shape.txt file
        outunit : str, optional
            Output units. ft or m, default is ft.

        Returns
        -------
        dfT : pandas.DataFrame
            Reformatted long-form DataFrame where each row represents a single sonar shot,
            including attributes such as depth, range, azimuth (cAzi), and tilt.
        """
        df = pd.read_csv(path, sep='\t', dtype=str)
        if 'H1 [m]' not in df.columns:
            raise Exception(path.replace('\\', '/'), 'Not a UMBRO shape format')

        df.columns = [x.replace(' [m]', '') for x in df.columns]  # abreviate that

        for col in list(df):
            df[col] = df[col].str.replace(',', '').astype(float)

        axi = np.linspace(0, 360, 17)[:-1]

        delta = [-11, -3.75, -7.5, 0.0, 3.75, 7.5, 11]  # 11.5 would be inbetween

        for i in range(1, 17):
            col = f'A{str(i)}'
            df[col] = axi[i - 1]

        dfT = pd.DataFrame()

        for i in range(1, 17):
            t = df[[x for x in df.columns if str(i) == x[1:]]].copy()
            t.columns = ['depth', 'r', 'cAzi']
            t['tilt'] = 0
            dfT = pd.concat([dfT, t], ignore_index=True)

            # Infill at equal radii not interplated
            for d in delta:
                td = t.copy()
                td.cAzi = td.cAzi + d
                dfT = pd.concat([dfT, td], ignore_index=True)

        ## Convert to ft?
        if outunit == 'ft':
            dfT['r'] = dfT['r'] * 3.28084
            dfT['depth'] = dfT['depth'] * 3.28084
    
        return dfT

    def infill_ubro_dft(dfT):
        """
        Infills depths within the dfT from top down.
        
        For each integer depth not present in the dataset, the function finds the 
        closest existing shallower depth and duplicates its data, updating only the 
        `depth` value to the missing level.

        Parameters
        ----------
        dfT : pandas.DataFrame
            A long-form DataFrame where each row represents a single sonar shot.
            Must contain columns `depth` and `cAzi`.

        Returns
        -------
        dfT : pandas.DataFrame
            Updated DataFrame with interpolated rows inserted at missing depths.
            Sorted in ascending order of depth.
        """
        '''
        Infills depths within the dfT from top down.
        '''
        depths = dfT.depth.values

        ## Infill per integer depth from top down
        aziz = dfT.cAzi.unique()
        for d in range(int(dfT.depth.min()), int(dfT.depth.max() + 1)):
            if d < dfT.depth.min():
                continue
            if d in dfT.depth.values:
                continue
            # Find the index of the point with the closest and less than values
            closest_index = np.argmin(np.abs(depths[depths <= d] - d))

            # Get the XYZ point with the closest X value
            closest_point = depths[closest_index]

            t = dfT[dfT.depth == closest_point].copy()
            t['depth'] = d
            dfT = pd.concat([dfT, t], ignore_index=True)

        dfT = dfT.sort_values('depth')

        return dfT