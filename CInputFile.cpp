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
#include "CInputFile.h"

//#include <QtCore>
#include <stdio.h>
#include <ogr_spatialref.h>
#include <iomanip>
#include <regex>
#include <cmath>

using namespace std;
extern "C"
{
    #include <jpeglib.h>
}

#define JPG_BLOCK_SIZE (TILESIZE*TILESIZE)

// this code is from the GDAL project
void printProgress(int current, int total)
{
    double dfComplete = double(current)/double(total);

    static int nLastTick = -1;
    int nThisTick = (int) (dfComplete * 40.0);

    nThisTick = MIN(40,MAX(0,nThisTick));
    
    // Have we started a new progress run? 
    if( nThisTick < nLastTick && nLastTick >= 39 )
        nLastTick = -1;
    
    if( nThisTick <= nLastTick )
        return;

    while( nThisTick > nLastTick ) {
        nLastTick++;
        if( nLastTick % 4 == 0 )
            fprintf( stdout, "%d", (nLastTick / 4) * 10 );
        else
            fprintf( stdout, "." );
    }

    if( nThisTick == 40 )
        fprintf( stdout, " - done.\n" );
    else
        fflush( stdout );
}



static std::vector<JOCTET> jpgbuf;
static void init_destination (j_compress_ptr cinfo) {
    jpgbuf.resize(JPG_BLOCK_SIZE);
    cinfo->dest->next_output_byte   = &jpgbuf[0];
    cinfo->dest->free_in_buffer     = jpgbuf.size();
}

static boolean empty_output_buffer (j_compress_ptr cinfo) { size_t oldsize = jpgbuf.size();
    jpgbuf.resize(oldsize + JPG_BLOCK_SIZE);
    cinfo->dest->next_output_byte   = &jpgbuf[oldsize];
    cinfo->dest->free_in_buffer     = jpgbuf.size() - oldsize;
    return true;
}
static void term_destination (j_compress_ptr cinfo) {
    jpgbuf.resize(jpgbuf.size() - cinfo->dest->free_in_buffer);
}


uint32_t CInputFile::nTilesTotal     = 0;
uint32_t CInputFile::nTilesProcessed = 0;

CInputFile::CInputFile(const string &filename, uint32_t tileSize, const string mapname, GDALRIOResampleAlg resampling_alg) {
    this->filename = filename;
    nTiles = 0;
    this->tileSize = tileSize;
    double adfGeoTransform[6]   = {0};
    char projstr[1024]          = {0};

    dataset = (GDALDataset*)GDALOpen(filename.c_str(),GA_ReadOnly);
    if(dataset == 0) {
        fprintf(stderr,"\nFailed to open %s\n", filename.c_str());
        exit(-1);
    }

    if(dataset->GetProjectionRef())
        strncpy(projstr,dataset->GetProjectionRef(),sizeof(projstr));
    char *p = projstr;
    oSRS.importFromWkt(&p);
    oSRS.exportToProj4(&p);
    pj = pj_init_plus(p);
    proj4 = p;
    CPLFree(p);

    if(pj == 0) {
        fprintf(stderr,"\nUnknown projection in file %s (%s)\n", filename.c_str(), projstr);
        exit(-1);
    }
    

    int32_t rasterBandCount = dataset->GetRasterCount();
    if(rasterBandCount == 1) {
        GDALRasterBand * pBand;
        pBand = dataset->GetRasterBand(1);


        if(pBand->GetColorInterpretation() ==  GCI_PaletteIndex ) {
            GDALColorTable * pct = pBand->GetColorTable();
            for(int i=0; i < pct->GetColorEntryCount(); ++i) {
                const GDALColorEntry& e = *pct->GetColorEntry(i);
                colortable[i] = ((e.c4 & 0x0ff) << 24) | ((e.c3 & 0x0ff) << 16) | ((e.c2 & 0x0ff) << 8) | (e.c1 & 0x0ff);
            }
        }
        else if(pBand->GetColorInterpretation() ==  GCI_GrayIndex ) {
            for(int i=0; i < 256; ++i)
                colortable[i] = i | (i << 8) | (i << 16) | 0xFF000000;
        }
        else {
            fprintf(stderr,"\nFile must be 8 bit palette or gray indexed.\n");
            exit(-1);
        }

        int success = 0;
        int idx = pBand->GetNoDataValue(&success);

        if(success)
            colortable[idx] &= 0x00FFFFFF;
    }

    dataset->GetGeoTransform( adfGeoTransform );

    width   = dataset->GetRasterXSize();
    height  = dataset->GetRasterYSize();

    if(pj_is_latlong(pj)) {
        xscale  = adfGeoTransform[1] * DEG_TO_RAD;
        yscale  = adfGeoTransform[5] * DEG_TO_RAD;
        xref1   = adfGeoTransform[0] * DEG_TO_RAD;
        yref1   = adfGeoTransform[3] * DEG_TO_RAD;
    } else {
        xscale  = adfGeoTransform[1];
        yscale  = adfGeoTransform[5];
        xref1   = adfGeoTransform[0];
        yref1   = adfGeoTransform[3];
    }
    xref2   = xref1 + width  * xscale;
    yref2   = yref1 + height * yscale;

    x1=xref1;
    y1=yref1;
    x2=xref2;
    y2=yref2;
    cropped_width=width;
    cropped_height=height;

    this->mapname = mapname;
    rasterio_args={1, resampling_alg, 0, 0, 0, 0,0,0,0};
}

CInputFile::~CInputFile()
{
}

void CInputFile::summarize() {
    printf("\n\n--- %s ---", filename.c_str());
    for(uint32_t i = 0; i < levels.size(); i++)
    {
        level_t& level = levels[i];
        printf("\nLevel%i:", i);
        printf("\nwidth/height:  %i/%i [pixel]", level.width, level.height);
        if(pj_is_latlong(pj))
        {
            printf("\nxscale/yscale: %1.6f/%1.6f [°/pixel]", level.xscale * RAD_TO_DEG, level.yscale * RAD_TO_DEG);
        }
        else
        {
            printf("\nxscale/yscale: %1.6f/%1.6f [m/pixel]", level.xscale, level.yscale);
        }

        printf("\nTiles X/Y:     %i/%i", level.xTiles, level.yTiles);
        printf("\n");
    }
}

void CInputFile::getRefP0(double& lon, double& lat) {
    projPJ wgs84 = pj_init_plus("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs");

    lon = x1;
    lat = y1;

    pj_transform(pj, wgs84, 1, 0, &lon, &lat, 0);
    pj_free(wgs84);

    lon *= RAD_TO_DEG;
    lat *= RAD_TO_DEG;
}

void CInputFile::getRefP1(double& lon, double& lat) {
    projPJ wgs84 = pj_init_plus("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs");

    lon = x2;
    lat = y1;

    pj_transform(pj, wgs84, 1, 0, &lon, &lat, 0);
    pj_free(wgs84);

    lon *= RAD_TO_DEG;
    lat *= RAD_TO_DEG;
}

void CInputFile::getRefP2(double& lon, double& lat) {
    projPJ wgs84 = pj_init_plus("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs");

    lon = x2;
    lat = y2;

    pj_transform(pj, wgs84, 1, 0, &lon, &lat, 0);
    pj_free(wgs84);

    lon *= RAD_TO_DEG;
    lat *= RAD_TO_DEG;
}

void CInputFile::getRefP3(double& lon, double& lat) {
    projPJ wgs84 = pj_init_plus("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs");

    lon = x1;
    lat = y2;

    pj_transform(pj, wgs84, 1, 0, &lon, &lat, 0);
    pj_free(wgs84);

    lon *= RAD_TO_DEG;
    lat *= RAD_TO_DEG;
}

void CInputFile::crop_to(CInputFile& base){
  OGRCoordinateTransformation * ct= OGRCreateCoordinateTransformation( &base.oSRS, &oSRS );
  x1=base.x1;
  y1=base.y1;
  x2=base.x2;
  y2=base.y2;
  if(!ct || !ct->Transform(1, &x1, &y1) || !ct->Transform(1, &x2, &y2)){
    fprintf(stderr, "Transformation failed.\n");
    exit(-1);
  }
  cropped_width = floor((x2 - x1) / xscale);
  cropped_height = floor((y2 - y1) / yscale);
}


uint32_t CInputFile::calcLevels(double scaleLimit, double& globXScale, double& globYScale) {
    nTiles     = 0;
    int32_t nLevels  = 1;

    if(scaleLimit) {
        double s = xscale;
        while(s*2.1 < scaleLimit) {
            s = s * 2;
            nLevels++;
        }
    } else {
        uint32_t w = cropped_width;
        while(w>>1 > TILESIZE) {
            w = w >> 1;
            nLevels++;
        }
    }

    for(int l = 0; l < nLevels; l++) {
        level_t level;

        level.xCorrectionScale = globXScale / xscale;
        level.yCorrectionScale = globYScale / yscale;
        level.xscale = (1 << l) * globXScale;
        level.yscale = (1 << l) * globYScale;
        level.width  = cropped_width *  xscale / level.xscale;
        level.height = cropped_height * yscale / level.yscale;
        level.xTiles = ceil(double(level.width) / tileSize);
        level.yTiles = ceil(double(level.height) / tileSize);
        level.number = 15-log2(level.xscale/5); // to do: what number does OruxMaps expect?

        levels.push_back(level);

        nTiles += level.xTiles * level.yTiles;

//        qDebug() << "level" << l << level.width << level.height << level.xTiles << level.yTiles << level.xscale << level.yscale;
    }

    globXScale = globXScale * (1 << nLevels);
    globYScale = globYScale * (1 << nLevels);

    nTilesTotal += nTiles;

    return nLevels;
}

void CInputFile::writeLevels(sqlite3 *db, std::ofstream& xmlstream, int quality, int subsampling) {
  for(uint32_t i = 0; i < levels.size(); i++) {
    writeLevelXml(xmlstream, i);
    writeLevel(db, i, quality, subsampling);
  }
}

void CInputFile::writeLevelXml(std::ofstream& xmlstream, int l){
    level_t& level = levels[l];

    double lon0, lat0;
    double lon1, lat1;
    double lon2, lat2;
    double lon3, lat3;
    getRefP0(lon0, lat0);
    getRefP1(lon1, lat1);
    getRefP2(lon2, lat2);
    getRefP3(lon3, lat3);

    // convert datum & projection from proj4 to oruxmaps strings
    // this needs more work
    string proj, datum, ellps;
    regex rgx;
    smatch match;
    rgx=".*\\+datum=(.*?) .*";
    if(regex_search(proj4, match, rgx))
      datum = match[1].str();
    rgx=".*\\+proj=(.*?) .*";
    if(regex_search(proj4, match, rgx))
      proj = match[1].str();

    if(proj == "tmerc"){
      proj = "Transverse Mercator";
      rgx=".*\\+lon_0=(.*?) .*";
      proj+=",";proj+=(regex_search(proj4, match, rgx))? match[1].str(): "0";
      rgx=".*\\+lat_0=(.*?) .*";
      proj+=",";proj+=(regex_search(proj4, match, rgx))? match[1].str(): "0";
      rgx=".*\\+k=(.*?) .*";
      proj+=",";proj+=(regex_search(proj4, match, rgx))? match[1].str(): "1";
      rgx=".*\\+x_0=(.*?) .*";
      proj+=",";proj+=(regex_search(proj4, match, rgx))? match[1].str(): "0";
      rgx=".*\\+y_0=(.*?) .*";
      proj+=",";proj+=(regex_search(proj4, match, rgx))? match[1].str(): "0";
    } else if(proj == "utm"){
      proj = "UTM,";
      rgx=".*\\+zone=(.*?) .*";
      if(regex_search(proj4, match, rgx)) proj+=match[1];
      else {
        fprintf(stderr, "UTM zone missing.\n");
        exit(-1);
      }
    } else {
      fprintf(stderr, "Unsupported projection: %s\n", proj.c_str());
      exit(-1);
    } 

    if(datum == "potsdam"){
      datum = "Potsdam Rauenberg DHDN:Poland@WGS 1984:Global Definition";
    } else {
      if(proj4.find("+ellps=GRS80") != string::npos)
        datum = "ETRS89:Europe@WGS 1984:Global Definition";
      else if(proj4.find("+ellps=bessel") != string::npos && proj4.find("+towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7"))
        datum = "Potsdam Rauenberg DHDN:Poland@WGS 1984:Global Definition";
      else {
        fprintf(stderr, "Unsupported datum: %s (%s)\n", datum.c_str(), proj4.c_str());
        exit(-1);
      }
    } 


    xmlstream << "<OruxTracker versionCode=\"2.1\">\n"
      <<"<MapCalibration layers=\"false\" layerLevel=\""<<level.number<<"\">\n"
      <<"<MapName><![CDATA["<<mapname<<" "<<level.number<<"]]></MapName>\n"
      <<"<MapChunks xMax=\""<< level.xTiles <<"\" yMax=\""<< level.yTiles <<"\" datum=\""<<datum<<"\" projection=\""<<proj<<"\" img_height=\"512\" img_width=\"512\" file_name=\""<<mapname<<" "<<level.number<<"\" />\n"
      <<"<MapDimensions height=\""<<level.height<<"\" width=\""<<level.width<<"\" />\n"
      <<"<MapBounds minLat=\""<<fixed<<setprecision(14)<<min(min(lat0, lat1),min(lat2,lat3))<<"\" maxLat=\""<<max(max(lat0, lat1),max(lat2,lat3))
        <<"\" minLon=\""<<min(min(lon0, lon1), min(lon2,lon3))<<"\" maxLon=\""<<max(max(lon0, lon1), max(lon2,lon3))<<"\" />\n"
      <<"<CalibrationPoints>\n"
      <<"<CalibrationPoint corner=\"TL\" lon=\""<<lon0<<"\" lat=\""<<lat0<<"\" />\n"
      <<"<CalibrationPoint corner=\"BR\" lon=\""<<lon2<<"\" lat=\""<<lat2<<"\" />\n"
      <<"<CalibrationPoint corner=\"TR\" lon=\""<<lon1<<"\" lat=\""<<lat1<<"\" />\n"
      <<"<CalibrationPoint corner=\"BL\" lon=\""<<lon3<<"\" lat=\""<<lat3<<"\" />\n"
      <<"</CalibrationPoints>\n"
      <<"</MapCalibration>\n"
      <<"</OruxTracker>\n";
}

void CInputFile::writeLevel(sqlite3 *db, int l, int quality, int subsampling) {
    level_t& level = levels[l];

    int32_t x0 = (x1-xref1)/xscale; // cropped
    int32_t y0 = (y1-yref1)/yscale;
    for(int y = 0; y < level.yTiles; y++) {
        for(int x = 0; x < level.xTiles; x++) {
            int32_t xoff = x0 + floor(x * (tileSize << l) * level.xCorrectionScale + 0.5);
            int32_t yoff = y0 + floor(y * (tileSize << l) * level.yCorrectionScale + 0.5);
            int32_t w1   = floor((tileSize << l) * level.xCorrectionScale + 0.5);
            int32_t h1   = floor((tileSize << l) * level.yCorrectionScale + 0.5);
            int32_t w2   = tileSize;
            int32_t h2   = tileSize;

            if(xoff + w1 > x0 + cropped_width) {
                w2 = w2 * (x0 + cropped_width - xoff) / w1;
                w1 = x0 + cropped_width - xoff;
            }
            if(yoff + h1 > y0 + cropped_height) {
                h2 = h2 * (y0 + cropped_height - yoff) / h1;
                h1 = y0 + cropped_height - yoff;
            }


     //cout <<l<<" "<<xoff<<" "<<yoff<<" "<<w1<<" "<<h1<<" "<<w2<<" "<<h2<<" "<<width<<" "<<height<<" "<<tileSize<<"\n";
            if(readTile(xoff, yoff, w1, h1, w2, h2, tileBuf32Bit)) {
                uint32_t size = compressTile(w2, h2, tileBuf32Bit, quality, subsampling);
                sqlite3_stmt *pStmt;
                const char sql[] = "INSERT INTO tiles VALUES(?, ?, ?, ?)";
                
                int rc = sqlite3_prepare(db, (char*)sql, -1, &pStmt, 0);
                if (rc != SQLITE_OK) {
                    fprintf(stderr, "Cannot prepare statement: %s\n", sqlite3_errmsg(db));
                    exit(-1);
                }    
                
                sqlite3_bind_int(pStmt, 1, x);
                sqlite3_bind_int(pStmt, 2, y);
                sqlite3_bind_int(pStmt, 3, level.number);
                sqlite3_bind_blob(pStmt, 4, (const void*)&jpgbuf[0], size, SQLITE_STATIC);
                
                rc = sqlite3_step(pStmt);
                
                if (rc != SQLITE_DONE){
                    fprintf(stderr, "execution failed: %s, %d, %d, %d, %d", sqlite3_errmsg(db),x,y,level.number, size);
                    exit(-1);
                }
                
                sqlite3_finalize(pStmt);
            } else {
                fprintf(stderr, "\nFailed to read tile from source\n");
                exit(-1);
            }

            nTilesProcessed++;
            printProgress(nTilesProcessed, nTilesTotal);
        }
    }
}

bool CInputFile::readTile(int32_t xoff, int32_t yoff, int32_t w1, int32_t h1, int32_t w2, int32_t h2, uint32_t *output) {
    int32_t rasterBandCount = dataset->GetRasterCount();

    //memset(output,-1, sizeof(uint32_t) * w2 * h2);

    if(rasterBandCount == 1) {
        GDALRasterBand * pBand;
        pBand = dataset->GetRasterBand(1);
        if(pBand->RasterIO(GF_Read, (int)xoff, (int)yoff, w1, h1, tileBuf08Bit, w2, h2, GDT_Byte,0,0,&rasterio_args) == CE_Failure)
            return false;

        for(int i = 0; i < (w2 * h2); i++)
            output[i] = colortable[(unsigned char)tileBuf08Bit[i]];
    } else {
        for(int b = 1; b <= rasterBandCount; ++b) {
            GDALRasterBand * pBand;
            pBand = dataset->GetRasterBand(b);

            uint32_t mask = ~(0x000000FF << (8*(b-1)));

            if(pBand->RasterIO(GF_Read,(int)xoff,(int)yoff, w1, h1, tileBuf08Bit,w2,h2,GDT_Byte,0,0,&rasterio_args) == CE_Failure)
                return false;

            for(int i = 0; i < (w2 * h2); i++) {
                uint32_t pixel = output[i];

                pixel &= mask;
                pixel |= tileBuf08Bit[i] << (8*(b-1));
                output[i] = pixel;
            }
        }
    }

    return true;
}

uint32_t CInputFile::compressTile(uint32_t xsize, uint32_t ysize, uint32_t * raw_image, int quality, int subsampling) {
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    JSAMPROW row_pointer[1];

    jpeg_destination_mgr destmgr    = {0};
    destmgr.init_destination        = init_destination;
    destmgr.empty_output_buffer     = empty_output_buffer;
    destmgr.term_destination        = term_destination;

    // convert from RGBA to RGB
    for(uint32_t r = 0; r < TILESIZE; r++)
        for(uint32_t c = 0; c < TILESIZE; c++) {
            uint32_t pixel;
            pixel = ((c>=xsize)||(r>=ysize))?0xffffff:raw_image[r * xsize + c];

            tileBuf24Bit[r * TILESIZE * 3 + c * 3]     =  pixel        & 0x0FF;
            tileBuf24Bit[r * TILESIZE * 3 + c * 3 + 1] = (pixel >>  8) & 0x0FF;
            tileBuf24Bit[r * TILESIZE * 3 + c * 3 + 2] = (pixel >> 16) & 0x0FF;
        }

    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress(&cinfo);

    cinfo.dest              = &destmgr;
    cinfo.image_width       = TILESIZE;//xsize;
    cinfo.image_height      = TILESIZE;//ysize;
    cinfo.input_components  = 3;
    cinfo.in_color_space    = JCS_RGB;

    jpeg_set_defaults( &cinfo );

    if (subsampling != -1) {
        switch (subsampling) {
        case 422:  // 2x1, 1x1, 1x1 (4:2:2) : Medium
            {
                cinfo.comp_info[0].h_samp_factor = 2;
                cinfo.comp_info[0].v_samp_factor = 1;
                cinfo.comp_info[1].h_samp_factor = 1;
                cinfo.comp_info[1].v_samp_factor = 1;
                cinfo.comp_info[2].h_samp_factor = 1;
                cinfo.comp_info[2].v_samp_factor = 1;
                break;
            }
        case 411:  // 2x2, 1x1, 1x1 (4:1:1) : High
            {
                cinfo.comp_info[0].h_samp_factor = 2;
                cinfo.comp_info[0].v_samp_factor = 2;
                cinfo.comp_info[1].h_samp_factor = 1;
                cinfo.comp_info[1].v_samp_factor = 1;
                cinfo.comp_info[2].h_samp_factor = 1;
                cinfo.comp_info[2].v_samp_factor = 1;
                break;
            }
        case 444:  // 1x1 1x1 1x1 (4:4:4) : None
            {
                cinfo.comp_info[0].h_samp_factor = 1;
                cinfo.comp_info[0].v_samp_factor = 1;
                cinfo.comp_info[1].h_samp_factor = 1;
                cinfo.comp_info[1].v_samp_factor = 1;
                cinfo.comp_info[2].h_samp_factor = 1;
                cinfo.comp_info[2].v_samp_factor = 1;
                break;
            }
        }
    }

    if (quality != -1)
        jpeg_set_quality( &cinfo, quality, TRUE );

    jpeg_start_compress( &cinfo, TRUE );

    while( cinfo.next_scanline < cinfo.image_height ) {
        row_pointer[0] = (JSAMPLE*)&tileBuf24Bit[ cinfo.next_scanline * cinfo.image_width *  cinfo.input_components];
        jpeg_write_scanlines( &cinfo, row_pointer, 1 );
    }
    
    jpeg_finish_compress( &cinfo );
    jpeg_destroy_compress( &cinfo );

    return jpgbuf.size();
}

