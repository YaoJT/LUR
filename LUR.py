## coding: utf-8

import arcpy
from sklearn import linear_model
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
import time

def extract(coor_list,raster_path):
##    no_value = arcpy.Describe(raster_path).noDataValue
    out_res = []
    array_f = arcpy.RasterToNumPyArray(raster_path)
    min_x = float(arcpy.GetRasterProperties_management(raster_path,'LEFT').getOutput(0))
    max_x = float(arcpy.GetRasterProperties_management(raster_path,'RIGHT').getOutput(0))
    min_y = float(arcpy.GetRasterProperties_management(raster_path,'bottom').getOutput(0))
    max_y = float(arcpy.GetRasterProperties_management(raster_path,'top').getOutput(0))
    cell_size_x =  float(arcpy.GetRasterProperties_management(raster_path,'cellsizex').getOutput(0))
    cell_size_y =  float(arcpy.GetRasterProperties_management(raster_path,'cellsizey').getOutput(0))

    for k in coor_list:
        point_x = k[0]
        point_y = k[1]
        if point_x >= min_x and point_x <= max_x and point_y >= min_y and point_y <= max_y:
            ncol = int((point_x-min_x)/cell_size_x)
            nrow = int((max_y-point_y)/cell_size_y)
            out_res.append(array_f[nrow,ncol])
        else:
            print "{0} is out of the extent({1}) of {2}".format(k,[min_x,max_x,min_y,max_y],raster_path)
            out_res.append(np.nan)
    return out_res

class LUR:
    """
    land_use map: the land use map or the path of land use map 
    factor_rasters: a list of other raster/raster path of the independent variables
    radius_range: a range of interested radiuse, such as range(1000,20000,1000)
    """
    def __init__(self,land_use_map,radius_range=range(1000,3000,1000),work_space='c:/LUR'):
        self.landuse_map = land_use_map
        self.radius_range = radius_range
        self.work_space = work_space
        arcpy.env.overwriteOutput=True
        if os.path.exists(self.work_space) == False:
            os.makedirs(self.work_space)

    def landmap_processing(self):
        out_path = self.work_space+'/land_map'
        if os.path.exists(out_path) == False:
            os.makedirs(out_path)
        land_map = []
        array = arcpy.RasterToNumPyArray(self.landuse_map)
        none_value = arcpy.Raster(self.landuse_map).noDataValue
        land_type = list(np.unique(array))
        land_type.remove(none_value)
        for i in land_type:
            map_path = out_path+'/land'+str(i)
            m = arcpy.sa.EqualTo(self.landuse_map,int(i))
            m.save(map_path)
            land_map.append(map_path)
        return land_map
    
    def normalize(self,rasters,method = 'min_max'):
        """
        method:
          Mean_std:
          min_max:
        """
        factor_map = []
        out_path = self.work_space+'/factor_map'
        if os.path.exists(out_path) == False:
            os.makedirs(out_path)
        for fac in rasters:
            map_path = os.path.join(out_path,os.path.split(fac)[1])
            array = arcpy.RasterToNumPyArray(fac)
            none_value = arcpy.Raster(fac).noDataValue
            value = []
            for i in array:
                for j in i:
                    if j != none_value:
                        value.append(j)
            array = np.array(value)
            if method.lower() == 'mean_std':
                mean_num,std_num = float(array.mean()),float(array.std())
                m = arcpy.sa.Minus(fac,mean_num)
                m = arcpy.sa.Divide(m,std_num)
            elif method.lower() == 'min_max':
                max_num,min_num = float(array.max()),float(array.min())
                m = arcpy.sa.Minus(fac,min_num)
                m = arcpy.sa.Divide(m,max_num-min_num)                
            m.save(map_path)
            factor_map.append(map_path)
        return factor_map

    def create_latlng(self):
        path = self.work_space+'/latlng_map'
        if os.path.exists(path) == False:
            os.makedirs(path)
            
        raster_path = self.landuse_map
        out_path = [path+'/lat',path+'/lng']
        array_f = arcpy.RasterToNumPyArray(raster_path)
        min_x = float(arcpy.GetRasterProperties_management(raster_path,'LEFT').getOutput(0))
        max_x = float(arcpy.GetRasterProperties_management(raster_path,'RIGHT').getOutput(0))
        min_y = float(arcpy.GetRasterProperties_management(raster_path,'bottom').getOutput(0))
        max_y = float(arcpy.GetRasterProperties_management(raster_path,'top').getOutput(0))
        cell_size_x =  float(arcpy.GetRasterProperties_management(raster_path,'cellsizex').getOutput(0))
        cell_size_y =  float(arcpy.GetRasterProperties_management(raster_path,'cellsizey').getOutput(0))
        num_row = int((max_y-min_y)/cell_size_y)+1
        num_col = int((max_x-min_x)/cell_size_x)+1
        lat_array,lng_array = np.ones((num_row,num_col)),np.ones((num_row,num_col))
        for i in range(num_row):
            for j in range(num_col):
                lat,lng = max_y-i*cell_size_y,min_x+j*cell_size_x
                lat_array[i][j],lng_array[i][j] = lat,lng
        lat_raster = arcpy.NumPyArrayToRaster(lat_array,arcpy.Point(min_x,lat),cell_size_x,cell_size_y)
        lng_raster = arcpy.NumPyArrayToRaster(lng_array,arcpy.Point(min_x,lat),cell_size_x,cell_size_y)
        lat_raster.save(out_path[0])
        lng_raster.save(out_path[1])
        return out_path
        
            

    def radius_identity(self,fit_array,maps,value_name = 'value',out_csv = True,fit_index = None):
        if fit_index == None:
            index = range(len(fit_array))
            print "fit_index has been automatically added since the length of fit_index is not equal to fit_array or not defined"              
        elif len(fit_index) == len(fit_array):
            index = fit_index
        else:
            index = range(len(fit_array))
            print "fit_index has been automatically added since the length of fit_index is not equal to fit_array or not defined"        
        column = ['X','Y',value_name]
        base_data = pd.DataFrame(fit_array,index,column)
        coor_list = []
        column_index = []
        for i in range(len(index)):
            coor_list.append((base_data.irow(i)['X'],base_data.irow(i)['Y']))
        out_path = self.work_space+'/radius_map'
        if os.path.exists(out_path) ==False:
            os.makedirs(out_path)
            
        for raster in maps:
            land = os.path.split(raster)[1].split('.')[0]
            for r in self.radius_range:
                neighbor = arcpy.sa.NbrCircle(r,'MAP')
                land_path = out_path+'/'+str(land)+'r_'+str(r)+'.tif'
                m = arcpy.sa.FocalStatistics(raster,neighbor,'Mean')
                m.save(land_path)
                value_list = extract(coor_list,land_path)
                base_data[land_path] = value_list
        for raster in maps:
            land = os.path.split(raster)[1].split('.')[0]
            corr = 0
            out_column = ''
            for i in self.radius_range:
                col = out_path+'/'+str(land)+'r_'+str(i)+'.tif'
                if abs(base_data[value_name].corr(base_data[col])) > corr+0.0000000001:
                    corr = abs(base_data[value_name].corr(base_data[col]))
                    out_column = col
            if out_column != '':
                column_index.append(out_column)
        if out_csv:
            base_data.to_csv(self.work_space+'/radius_identity.csv')
            base_data[column_index].to_csv(self.work_space+'/radius_identity_result.csv')
        return column_index
                
    def extract_values(self,fit_array,factor_maps,value_name = 'value',out_csv = True,fit_index = None):
        if fit_index == None:
            index = range(len(fit_array))
            print "fit_index has been automatically added since the length of fit_index is not equal to fit_array or not defined"              
        elif len(fit_index) == len(fit_array):
            index = fit_index
        else:
            index = range(len(fit_array))
            print "fit_index has been automatically added since the length of fit_index is not equal to fit_array or not defined" 
        
        column = ['X','Y',value_name]
        base_data = pd.DataFrame(fit_array,index,column)
        coor_list = []
        for i in range(len(index)):
            coor_list.append((base_data.irow(i)['X'],base_data.irow(i)['Y']))        

## extract value to df
        for name in factor_maps:
            value_list = extract(coor_list,name)
            base_data[name] = value_list
            base_data.to_csv(os.path.join(self.work_space,'extract_values_'+value_name+'.csv'))
        return base_data       

    def factor_identity(self,base_data,out_csv = True,min_corr = 0.3,fix_rasters = []):
## calculate correlation coefficient
        dff = base_data.icol(range(2,len(base_data.columns))).corr()
## condition: corr > min_corr
        column_index = []
        for i in range(len(dff.columns)):
            if (abs(dff.icol(0).irow(i)) > min_corr) or (dff.index[i] in fix_rasters):
                column_index.append(dff.index[i])
        out_df = base_data[column_index]
##        print column_index
## writing
        if out_csv:
            dff.to_csv(os.path.join(self.work_space,'factor_identity_corr.csv'))
            out_df.to_csv(os.path.join(self.work_space,'factor_identity_result.csv'))
        return out_df
            
        
    def fit(self,fit_df,component_variable=0.85,fit_intercept = True,LOOCV = True):
        X_num = len(fit_df.columns)
        num = len(fit_df.index)        
        index = fit_df.index
        pca = PCA(n_components = component_variable)
        clf = linear_model.LinearRegression()
        clf.fit_intercept = fit_intercept
        new_X = pca.fit_transform(fit_df.icol(range(1,X_num)))
        data_Y = fit_df.icol(0)
        clf.fit(new_X,data_Y)

        pca_out,clf_out = {},{}
        help_x = np.eye(X_num-1)
        for i in range(1,X_num):
            pca_out[fit_df.columns[i]]= pca.transform(help_x[i-1])[0]
        if fit_intercept:
            intercept = clf.predict([0]*len(new_X[0]))
            help_x = np.eye(len(new_X[0]))
            for i in range(len(new_X[0])):
                clf_out[i] = clf.predict(help_x[i])-intercept
            clf_out['intercept'] = intercept
        else:
            help_x = np.eye(num_component)
            for i in range(num_component):
                clf_out[i] = clf.predict(help_x[i])
            clf_out['intercept'] = 0
## add statistical indexes for pca and clf (unfinished)
        out_file = open(self.work_space+'/pca_result.csv','w+')
        out_file.write('pca\n')


        out_file.write('regression\n')


        out_file.write('statistical values\n')


        out_file.close()
## next have been finished        
        if LOOCV:
            cv_out = 0
            out_f = os.path.join(self.work_space,'LOOCV.csv')
            out_f = open(out_f,'w+')
            out_f.write('time,{0},component_variable,{1}\n'.format(time.ctime(),component_variable))
            out_f.write('index,real,predict\n')
            for i in range(num):
                out_f.write(index[i]+',')
                real_X = new_X[i]
                real_y = data_Y[i]
                fit_X = [new_X[x] for x in range(num) if x != i]
                fit_y = [data_Y[x] for x in range(num) if x != i]
                clf = linear_model.LinearRegression()
                clf.fit_intercept = fit_intercept
                clf.fit(fit_X,fit_y)
                pre_y = clf.predict(real_X)
                out_f.write(str(real_y)+',')
                out_f.write(str(pre_y)+'\n')
                cv_out += (pre_y-real_y)**2
            out_f.close()
            cv_out = cv_out**0.5/num
            return pca_out,clf_out,cv_out
        else:
            return pca_out,clf_out

    def optimize(self,fit_df,component_variable=0.85,min_corr=0.3,fit_intercept = True, step_com = 0.05,step_corr = 0.05,fix_rasters = []):
        fix = fix_rasters    
        min_mse = 999999999999
        optm_com = component_variable
        optm_corr = min_corr
        out_file = open(self.work_space+'/optimize.csv','w+')
        out_file.write('corr,component,MSE\n')
        for corr in np.arange(min_corr,1,step_corr):
            fit_data = self.factor_identity(fit_df,min_corr= corr,fix_rasters=fix)
            if len(fit_data.columns)>2:
                for component in np.arange(component_variable,1,step_com):
                    mse = self.fit(fit_data,component,fit_intercept)[2]
                    out_file.write('{0},{1},{2}\n'.format(corr,component,mse))
                    if mse < min_mse-0.00000001:                    
                        min_mse = mse
                        print 'optimizing: coor: {0}, component: {1}, mse: {2}'.format(corr,component,min_mse)
                        optm_com,optm_corr = component,corr
        out_file.write('result\n')
        out_file.write('{0},{1},{2}\n'.format(optm_corr,optm_com,min_mse))
        fit_data = self.factor_identity(fit_df,min_corr= optm_corr,fix_rasters=fix)
        pca,clf = self.fit(fit_data,component_variable=optm_com,fit_intercept=fit_intercept)[:2]
        out_file.close()
        return pca,clf,min_mse
        

    def mapping(self,clf,pca,out_name = 'mapping_result',mask = None):
        if mask != None:
            arcpy.env.mask = mask
            arcpy.env.extent = mask
        path = self.work_space+'/components'
        if os.path.exists(path) == False:
            os.makedirs(path)
        components = []
        out_raster = clf['intercept']
        for i in range(len(clf)-1):
            out_path = os.path.join(path,'component_'+str(i))
            res = 0
            for factor in pca:
                res = arcpy.sa.Plus(res,arcpy.sa.Times(factor,pca[factor][i]))
            out_raster = arcpy.sa.Plus(out_raster,arcpy.sa.Times(res,clf[i]))
            res.save(out_path)
        
            
        out_raster.save(os.path.join(self.work_space,out_name))
        return os.path.join(self.work_space,out_name)


if __name__ == '__main__':
    optimize = False
    path = 'g:/GHG_Landuse/python_module/text_data/data.txt'
    fi = open(path).readlines()
    data = [f.replace('\n','').split() for f in fi]
    column = data[0]
    index = [da[0] for da in data[1:]]
    for p in range(7,len(data[0])):
        value_list = [da[1:3]+[da[p]] for da in data[1:]]    
##        value_list = [da[1:] for da in data]
        value_name = column[p]
        for i in range(len(value_list)):
            for j in range(len(value_list[i])):
                try:
                    value_list[i][j] = float(value_list[i][j])
                except:
                    None
        print "start processing for "+value_name
       
        arcpy.CheckOutExtension('spatial')
        land_map = 'g:/GHG_Landuse/python_module/text_data/lucc2015_f.tif'
        factors = ['g:/GHG_Landuse/python_module/text_data/dem.tif','g:/GHG_Landuse/python_module/text_data/slope.tif',
                   'g:/GHG_Landuse/python_module/data/dis_road','g:/GHG_Landuse/python_module/data/dis_water']
        lur = LUR(land_map,work_space = 'h:/lur2_'+value_name,radius_range=range(1000,31000,2000))
## procesing land maps and creater latlng maps
        land_maps = lur.landmap_processing()
        other_maps = lur.create_latlng()
## normalization
        factor_maps1 = lur.normalize(factors)
        factor_maps2 = lur.normalize(other_maps)
## radius identity
        radius_maps = lur.radius_identity(value_list,land_maps+factor_maps1,fit_index=index,value_name = value_name)
## get factor maps
        factor_maps = radius_maps + factor_maps2
## extract values
        base_data = lur.extract_values(value_list,factor_maps,fit_index=index,value_name = value_name)
## 
        if optimize:
            pca,clf = lur.optimize(base_data,component_variable=0.8,min_corr = 0.3)[:2]
        else:
            fit_data = lur.factor_identity(base_data,min_corr = 0.35)
            pca,clf = lur.fit(fit_data,component_variable=0.85)[:2]
            
        lur.mapping(clf,pca,value_name,mask = factors[0])
        arcpy.CheckInExtension('spatial')
        break
