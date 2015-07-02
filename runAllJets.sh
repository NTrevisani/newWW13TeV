
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi

LUMINOSITY=5.

NJETS=$1

CHANNELS="All SF OF EE MuE EMu MuMu "

PROOFMODE="Sequential"

SAMESIGN="OS" 

SAMPLES="
WW50                \
WJets50             \
TTbar50             \
"
#QCD                \ 
#Top                \
#TTJets             \
#VBF                \
#WW                 \

#rm -rf rootfiles/${NJETS}jet

mkdir rootFiles

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	mkdir rootFiles/
	mkdir rootFiles/${CHANNEL}	
	root -l -b -q "RunMuonAnalyzer.C(\"$SAMPLE\",\"$CHANNEL\",\"$SAMESIGN\",\"$PROOFMODE\",$LUMINOSITY)"
  
    done

done