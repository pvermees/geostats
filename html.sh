rm -rf html/*
cp -r figures html
cp -r latex html
cd html/latex
pdflatex geostats.tex
pdflatex geostats.tex
make4ht -d .. geostats.tex
cd ../..
