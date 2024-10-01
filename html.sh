rm -rf html/*
cp -r figures html
cp -r latex html
cd html/latex
rm *.log
rm *.aux
rm *.out
rm *.toc
pdflatex ./geostats.tex
pdflatex ./geostats.tex
pdflatex ./geostats.tex
make4ht ./geostats.tex
