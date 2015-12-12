/**********************************************************************************************
    Copyright (C) 2012 Oliver Eichler oliver.eichler@gmx.de
                  2015 Adrian Hänsler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************************/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <gdal_priv.h>
#include <proj_api.h>
#include <ogr_spatialref.h>
#include <sqlite3.h>

#define TILESIZE 512

extern void printProgress(int current, int total);

class CInputFile
{
    public:
        CInputFile(const std::string& filename, uint32_t tileSize, const std::string mapname, GDALRIOResampleAlg resampling_alg);
        virtual ~CInputFile();

        void summarize();

        double getXScale(){return xscale;}
        double getYScale(){return yscale;}

        void getRefP0(double& lon, double& lat);
        void getRefP1(double& lon, double& lat);
        void getRefP2(double& lon, double& lat);
        void getRefP3(double& lon, double& lat);

        int32_t getWidth(){return width;}
        int32_t getHeight(){return height;}
        std::string getProjection(){return compeProj;}
        std::string getDatum(){return compeDatum;}

        void crop_to(CInputFile& base);
        uint32_t calcLevels(double scaleLimit, double& globXScale, double& globYScale);

        void writeLevels(sqlite3 *db, std::ofstream& xmlstream, int quality, int subsampling);

        static uint32_t getTilesTotal(){return nTilesTotal;}

        std::string mapname;
        OGRSpatialReference oSRS;
        double xscale;
        double yscale;
        double xref1;
        double yref1;
        double xref2;
        double yref2;
        double x1, y1;
        double x2, y2;
        uint32_t cropped_width;
        uint32_t cropped_height;

    private:
        std::string proj4;
        void writeLevelXml(std::ofstream& xmlstream, int level);
        void writeLevel(sqlite3 *db, int level, int quality, int subsampling);

        bool readTile(int32_t xoff, int32_t yoff, int32_t w1, int32_t h1, int32_t w2, int32_t h2, uint32_t *output);
        uint32_t compressTile(uint32_t xsize, uint32_t ysize, uint32_t * raw_image, int quality, int subsampling);

        std::string filename;
        uint32_t tileSize;

        struct level_t
        {
            level_t(): number(0), width(0), height(0), xTiles(0), yTiles(0), xscale(0), yscale(0), xCorrectionScale(0), yCorrectionScale(0){}
            int32_t number;
            int32_t width;
            int32_t height;
            int32_t xTiles;
            int32_t yTiles;
            //QVector<quint64> offsetJpegs;

            double xscale;
            double yscale;

            double xCorrectionScale;
            double yCorrectionScale;
        };

        projPJ  pj;
        GDALDataset * dataset;
        std::string compeProj;
        std::string compeDatum;

        int32_t width;
        int32_t height;


        uint32_t nTiles;
        static uint32_t nTilesTotal;
        static uint32_t nTilesProcessed;
        std::vector<level_t> levels;
        char tileBuf08Bit[TILESIZE * TILESIZE];
        char tileBuf24Bit[TILESIZE * TILESIZE * 3];
        uint32_t tileBuf32Bit[512*512];
        uint32_t colortable[256];
        GDALRasterIOExtraArg rasterio_args;
};

