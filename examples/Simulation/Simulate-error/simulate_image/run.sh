./0mtz_to_P1hkl.com P43212.mtz

./simulate-errors-p43212

indexamajig -i pat.lst -g simple.geom --indexing=none --peaks=zaef \
--min-gradient=90  --min-snr=0.2 --threshold=10 \
-o peaks.stream

#./awk.sh
