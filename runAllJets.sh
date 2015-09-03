
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi

LUMINOSITY=0.04008

NJETS=$1

CHANNELS="MuMu" 
#MuMu All SF EE MuE EMu MuMu "

PROOFMODE="Cluster"

SAMESIGN="OS" 

MUONIDS="MediumIDTighterIP"
#"MediumID TightID TightIDTighterIP MediumIDTighterIP"

SAMPLES="
Top50                \
Data2015_50          \
WW50                 \
VV50                 \
WJets50              \
DY50                 \
"
#WW25                 \
#WJets25              \
#HWW25                \
#ZZ25                 \
#singleTop25          \
#TW25                 \
#"

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