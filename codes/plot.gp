# Live 2D colormap of FT2.dat

# Use an interactive terminal (wxt or qt)
set terminal wxt size 800,600 enhanced font 'Arial,14'
# OR: set terminal qt

# Title and labels
set title "2D Colormap of FT2.dat"
set xlabel "j"
set ylabel "i"

# Color map settings
set pm3d map                # top-down view
set palette rgb 33,13,10    # smooth palette
set colorbox                # show color scale

# Optional: set axis and color scale ranges
# set xrange [0:50]
# set yrange [40:60]
#set cbrange [0.5:0.7]

# Plot
splot '../data/rho2Deps1.000000ew1.000000dx0.025000rhob0.597846.dat' using 2:1:3 with pm3d notitle
pause -1
