
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi

LUMINOSITY=0.04008

NJETS=$1

CHANNELS="SF OF MuMu EE"
#"MuMu All SF EE MuE EMu MuMu "

PROOFMODE="Lite"

SAMESIGN="OS" 

MUONIDS="MediumIDTighterIP"
#"MediumIDTighterIP MediumID TightID TightIDTighterIP"

SAMPLES="
HWW50                \
Data201550           \
Top50                \
WW50                 \
VV50                 \
WJets50              \
DY25                 \
DY50                 \
TTJets50             \
"

mkdir rootFiles

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 

	for MUONID in $MUONIDS; do 
	
	    mkdir rootFiles/
	    mkdir rootFiles/${CHANNEL}	
	    mkdir rootFiles/${CHANNEL}/${MUONID}	
	    root -l -b -q "RunMuonAnalyzer.C(\"$SAMPLE\",\"$CHANNEL\",\"$SAMESIGN\",\"$PROOFMODE\",$LUMINOSITY,\"$MUONID\")"
	    
	    endproof -f *
	
	    resetpaf
	
	done
	
    done

done

#hadd -f rootFiles/${CHANNEL}/${MUONID}/Top.root rootFiles/${CHANNEL}/${MUONID}/TTJets50.root rootFiles/${CHANNEL}/${MUONID}/Top50.root
hadd -f rootFiles/${CHANNEL}/${MUONID}/DY.root rootFiles/${CHANNEL}/${MUONID}/50ns/DY.root rootFiles/${CHANNEL}/${MUONID}/25ns/DY.root
mv rootFiles/${CHANNEL}/${MUONID}/DY.root rootFiles/${CHANNEL}/${MUONID}/50ns/DY.root