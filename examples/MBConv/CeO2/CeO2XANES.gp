shift=0.46
set xlab 'Energy (eV)'
set ylab 'Absorption (arb. units)'
p 'CeO2M45.exp' u ($1):($2/5.0+0.2) w lp tit 'Expt.', \
'mbxanes.dat' u ($1-shift):($2+0.1) tit 'Many-body', \
'' u ($1-shift):($3/2.5) tit 'Single particle'
pause -1
set out 'CeO2.ps'
@psterm
rep
