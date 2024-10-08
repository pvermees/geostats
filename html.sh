rm -rf html/*
cp -r figures html
cp -r latex html
cd html/latex
rm ./*.log
rm ./*.aux
rm ./*.out
rm ./*.toc
find . -name "*.tex" -exec sed -i \
     's/\\includegraphics\[width=\\textwidth\]/\\includegraphics\[\]/g' {} +
find . -name "*.tex" -exec sed -i \
     's/\\includegraphics\[width=\\linewidth\]/\\includegraphics\[\]/g' {} +
find . -name "*.tex" -exec sed -i 's/\\boxed/\\fbox/g' {} +
pdflatex ./geostats.tex
pdflatex ./geostats.tex
pdflatex ./geostats.tex
make4ht ./geostats.tex
sed -i 's/<\/head>/<style>\nbody\{width:800px;margin-left:auto;margin-right:auto;\}\nimg\{max-height:200px; max-width:800px\;\}\n<\/style>\n<\/head>/g' ./geostats.html
sed -i 's/\.cmr\-6{/\.cmr\-6{padding\-right:1em;/g' ./geostats.css
mv ./geostats.html ./index.html
