#include <iostream>
#include <fstream>

#include "model/model.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MEpoch.h>

#include <StationResponse/LofarMetaDataUtil.h>
#include <StationResponse/ITRFConverter.h>

#include <aocommon/banddata.h>
#include <aocommon/matrix2x2.h>

using aocommon::MC2x2;

void dirToITRF(LOFAR::StationResponse::ITRFConverter& converter, const casacore::MDirection& dir, LOFAR::StationResponse::vector3r_t& itrf)
{
	casacore::MDirection itrfDir = converter.toDirection(dir);
	casacore::Vector<double> itrfVal = itrfDir.getValue().getValue();
	itrf[0] = itrfVal[0];
	itrf[1] = itrfVal[1];
	itrf[2] = itrfVal[2];
}

void sourceResponse(const std::string& msFilename, const ModelComponent& source, const std::string& filename)
{
  casacore::MeasurementSet ms(msFilename);
  
  aocommon::BandData band(ms.spectralWindow());
  double subbandFrequency = band.CentreFrequency();
  
  casacore::MSField fieldTable(ms.field());
  if(fieldTable.nrow() != 1)
    throw std::runtime_error("Set has multiple fields");
  casacore::ScalarMeasColumn<casacore::MDirection> delayDirColumn(fieldTable, casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));
  casacore::MDirection delayDir = delayDirColumn(0);
  
  casacore::MDirection tileBeamDir;
  if(fieldTable.tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
    casacore::ROArrayMeasColumn<casacore::MDirection> tileBeamDirColumn(fieldTable, "LOFAR_TILE_BEAM_DIR");
    tileBeamDir = *(tileBeamDirColumn(0).data());
  } else {
    tileBeamDir = delayDir;
  }
  std::vector<LOFAR::StationResponse::Station::Ptr> stations(ms.antenna().nrow());
  readStations(ms, stations.begin());
  
  casacore::MEpoch::ScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
  double startTime = timeColumn(0).getValue().get()*86400.0;
  double prevTime = startTime - 1.0;
  
  std::ofstream file(filename);
  for(size_t row=0; row!=ms.nrow(); ++row)
  {
    casacore::MEpoch time = timeColumn(row);
    double timeAsDouble = time.getValue().get()*86400.0;
    if(timeAsDouble != prevTime)
    {
      prevTime = timeAsDouble;
      LOFAR::StationResponse::ITRFConverter itrfConverter(timeAsDouble);
      
      LOFAR::StationResponse::vector3r_t station0, tile0;
      dirToITRF(itrfConverter, delayDir, station0);
      dirToITRF(itrfConverter, tileBeamDir, tile0);
      
      LOFAR::StationResponse::vector3r_t itrfDirection;
      static const casacore::Unit radUnit("rad");
      casacore::MDirection sourceDir(casacore::MVDirection(
        casacore::Quantity(source.PosRA(), radUnit),
        casacore::Quantity(source.PosDec(),radUnit)),
        casacore::MDirection::J2000);
      
      double stokesI = source.SED().FluxAtFrequency(subbandFrequency, aocommon::Polarization::StokesI);

      dirToITRF(itrfConverter, sourceDir, itrfDirection);
      MC2x2 response = MC2x2::Zero();
      size_t count = 0;
      double maxEigenValue = 0.0;
      for(size_t station=0; station!=stations.size(); ++station)
      {
        LOFAR::StationResponse::matrix22c_t gainMatrix =
          stations[station]->response(timeAsDouble, subbandFrequency, itrfDirection, subbandFrequency, station0, tile0);
        
        MC2x2 stationResponse( gainMatrix[0][0], gainMatrix[0][1], gainMatrix[1][0], gainMatrix[1][1] );
        response += stationResponse;
        ++count;
        std::complex<double> e1,e2;
        (stationResponse * stokesI).EigenValues(e1, e2);
        maxEigenValue = std::max(maxEigenValue, std::max(std::abs(e1), std::abs(e2)));
      }
      response = response * (1.0 / count);
      std::complex<double> e1, e2;
      (response * stokesI).EigenValues(e1, e2);
      file << ((timeAsDouble-startTime)/3600.0) << '\t' << maxEigenValue << '\t' << std::max(std::abs(e1), std::abs(e2)) << '\n';
    }
  }
}

std::ofstream header(const std::string& name)
{
  std::ofstream responsePlt(name + ".plt");
  responsePlt <<
    "set terminal pdfcairo enhanced color font 'Times,16'\n"
    "#set logscale xy\n"
    "#set yrange [0.1:]\n"
    "set output \"" << name << ".pdf\"\n"
    "set key above\n"
    "set xlabel \"Time (h)\"\n"
    "set ylabel \"Apparent flux (Jy)\"\n"
    "plot \\\n";
  return responsePlt;
}

int main(int argc, char* argv[])
{
  if(argc < 3)
    std::cout << "Syntax: sourceresponse <ms> <model>\n";
  
  std::ofstream maxPlt(header("response-max"));
  std::ofstream avgPlt(header("response-avg"));

  Model model(argv[2]);
  bool isFirst = true;
  for(const ModelSource& s : model)
  {
    for(size_t i=0; i!=s.ComponentCount(); ++i)
    {
      const ModelComponent& c = s.Component(i);
      std::string name = (s.ComponentCount()!=1) ? s.Name() + "_" + std::to_string(i) : s.Name();
      std::cout << "Calculating " << name << "...\n";
      sourceResponse(argv[1], c, name + ".txt");
      if(!isFirst) {
        maxPlt << ",\\\n";
        avgPlt << ",\\\n";
      }
      else {
        isFirst = false;
      }
      maxPlt << '"' << name << ".txt\" using 1:2 with lines title '" << name << "' lw 2";
      avgPlt << '"' << name << ".txt\" using 1:3 with lines title '" << name << "' lw 2";
   }
  }
  maxPlt << "\n";
  avgPlt << "\n";
}

