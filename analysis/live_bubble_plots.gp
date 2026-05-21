# ==============================================================================
# LIVE SIMULATION MONITORING SCRIPT USING GNUPLOT
# ==============================================================================
# 
# WHAT THIS SCRIPT DOES:
# This script creates a real-time, live-updating graphical dashboard for the
# simulation data. It targets the file "bubbleV.dat" and extracts 3 different
# physical metrics, plotting them simultaneously against Time (Column 1).
#
# EXPECTED DATA FORMAT (bubbleV.dat):
# The script expects a file with at least 4 columns:
# Column 1: Time (X-Axis for all plots)
# Column 2: totalHe         Column 3: totalV          Column 4: bubbleDensity
#
# HOW TO RUN IT: 
# gnuplot live_bubble_plots.gp
#
# HOW TO INTERACT WITH THE PLOT:
# - This is a hands-off, automated viewer. It will scale itself and refresh
#   automatically every second as new data is written by the simulation.
# - Close: Press 'q' on your keyboard or Ctrl+C in the terminal to exit.
#
# PREREQUISITES:
# Requires 'gnuplot-qt' or 'gnuplot-x11' for graphical display windows.
# ==============================================================================

set terminal qt size 800,1000  # Sets a tall window for stacked plots

# Start an infinite loop using the "while" syntax
while (1) {
    set multiplot layout 3,1 title "Live Simulation Data" font ",14"

    # Common settings for all plots
    set xlabel "Time"
    set grid

    # Plot 1: totalHe
    set ylabel "totalHe"
    plot "bubbleV.dat" using 1:2 with lines lc rgb "blue" notitle

    # Plot 2: totalV
    set ylabel "totalV"
    plot "bubbleV.dat" using 1:3 with lines lc rgb "red" notitle

    # Plot 3: bubbleDensity
    set ylabel "bubbleDensity"
    plot "bubbleV.dat" using 1:4 with lines lc rgb "green" notitle

    unset multiplot
    
    # Pause for 1 second before looping and replotting
    pause 1
}
