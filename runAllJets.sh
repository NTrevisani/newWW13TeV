
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi

LUMINOSITY=5.

NJETS=$1

CHANNELS="OF"
#All SF EE MuE EMu MuMu "

PROOFMODE="Cluster"

SAMESIGN="OS" 

MUONIDS="MediumIDTighterIP"
#MediumID TightID TightIDTighterIP"

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

	for MUONID in $MUONIDS; do 
	
	    mkdir rootFiles/
	    mkdir rootFiles/${CHANNEL}	
	    mkdir rootFiles/${CHANNEL}/${MUONID}	
	    root -l -b -q "RunMuonAnalyzer.C(\"$SAMPLE\",\"$CHANNEL\",\"$SAMESIGN\",\"$PROOFMODE\",$LUMINOSITY,\"$MUONID\")"
	    
	done
	
    done

done