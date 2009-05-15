#! /bin/sh
# the next line restarts using wish \
exec wish $0 ${1+"$@"}

package require BLT

########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################

#
# Graphical User Interface for the "sizes" program
#  edited September 1999
#

######################################################################
#
#   doSasPlot
#      plots the SAS fit from an analysis
#
proc doSasPlot {} {
	global params plotnum

	#
	# define the local variables
	#
	set ext fit
	set win plot_$plotnum
	set dataFile $params(wdName)/$params(project).$ext
	set graphWindow .$win:$ext
	set graph $graphWindow.graph

	#
	# read the data file
	#
	set Q   {}
	set I   {}
	set FIT {}
	set ESD {}
	set Z   {}
	set f [open $dataFile r]
	set line [ gets $f ];			# title line
		set title [string range $line 11 end]
	gets $f;			# column headings
	while {![eof $f]} {
		set line [ gets $f ]
		if {[string length $line] == 0} break;
		scan $line "%s%s%s%s%s" theQ theI theESD theFIT theZ
		if {$theI > 0} {
			lappend Q   $theQ
			lappend I   $theI
			lappend ESD $theESD
			lappend FIT $theFIT
			lappend Z   $theZ
		}
	}
	close $f

	#
	# use BLT utilities for Tcl/Tk
	#
	if {[winfo exists $graphWindow]} {destroy $graphWindow}
	toplevel $graphWindow
	blt::graph $graph
	$graph configure -title $title -font $params(font)
	$graph element create $params(sasFile) -xdata $Q -ydata $I \
	    -symbol circle -linewidth 0 -pixels 0.5m
	$graph element create fit -xdata $Q -ydata $FIT \
	    -symbol none -pixels 1m -color red
	$graph legend configure -position top -hide no -font $params(font)
	$graph crosshairs configure -hide no
	$graph axis configure x -title "Q, 1/Angstroms"
	$graph axis configure y -title "intensity, 1/cm"
	$graph crosshairs configure -linewidth 1 -hide no
	$graph axis configure x -logscale true
	$graph axis configure y -logscale true

	pack append $graphWindow \
		$graph { fill expand padx 20 pady 20 }

	wm min $graphWindow 10 10

	$graph postscript output $params(project)-$ext.ps \
	        -colormode color -width 5i -height 5i

        Blt_ZoomStack    $graph
        Blt_Crosshairs   $graph
        Blt_ActiveLegend $graph
        Blt_ClosestPoint $graph

	incr plotnum
}

######################################################################
#
#   doResPlot
#      plots the residuals from an analysis
#
proc doResPlot {} {
	global params plotnum

	#
	# define the local variables
	#
	set ext fit
	set win plot_$plotnum
	set dataFile $params(wdName)/$params(project).$ext
	set graphWindow .$win:$ext
	set graph $graphWindow.graph

	#
	# read the data file
	#
	set Q {}
	set Z {}
	set f [open $dataFile r]
	set line [ gets $f ];			# title line
		set title [string range $line 11 end]
	gets $f;			# column headings
	while {![eof $f]} {
		set line [ gets $f ]
		if {[string length $line] == 0} break;
		scan $line "%s%*s%*s%*s%s" theQ theZ
		lappend Q $theQ
		lappend Z $theZ
	}
	close $f

	#
	# use BLT utilities for Tcl/Tk
	#
	if {[winfo exists $graphWindow]} {destroy $graphWindow}
	toplevel $graphWindow
	blt::graph $graph
	$graph configure -title $title -font $params(font)
	$graph element create line1 -xdata $Q -ydata $Z \
	    -symbol circle -linewidth 0 -pixels 1m
	$graph legend configure -hide yes
	$graph crosshairs configure -hide no
	$graph axis configure x -title "Q, 1/Angstroms"
	$graph axis configure y -title "standardized residual"
	#$graph axis configure x -logscale true

	pack append $graphWindow \
		$graph { fill expand padx 20 pady 20 }

	wm min $graphWindow 10 10

	$graph postscript output $params(project)-resid.ps \
	        -colormode color -width 5i -height 5i

        Blt_ZoomStack    $graph
        Blt_Crosshairs   $graph
        Blt_ActiveLegend $graph
        Blt_ClosestPoint $graph
	incr plotnum
}

######################################################################
#
#   doSizPlot
#      plots the size distribution from an analysis
#
proc doSizPlot {} {
	global params plotnum
	if {$params(distType) == 0} {set ext "N-D"; set lbl "N(D)"}
	if {$params(distType) == 1} {set ext "f-D"; set lbl "f(D), 1/A"}
	if {$params(distType) == 2} {set ext "i-D"; set lbl "i(D)"}

	#
	# define the local variables
	#
	set win plot_$plotnum
	set dataFile $params(wdName)/$params(project).$ext
	set graphWindow .$win:$ext
	set graph $graphWindow.graph

	#
	# read the data file
	#
	set D   {}
	set QTY {}
	set f [open $dataFile r]
	set line [ gets $f ];			# title line
		set title [string range $line 11 end]
	gets $f;			# column headings
	while {![eof $f]} {
		set line [ gets $f ]
		if {[string length $line] == 0} break;
		scan $line "%s%s" theD theQTY
		lappend D   $theD
		lappend QTY $theQTY
	}
	close $f

	#
	# use BLT utilities for Tcl/Tk
	#
	if {[winfo exists $graphWindow]} {destroy $graphWindow}
	toplevel $graphWindow
	blt::graph $graph
	$graph configure -title $title -font $params(font)
	$graph element create line1 -xdata $D -ydata $QTY \
	    -symbol circle -linewidth 1 -pixels 1m
	$graph legend configure -hide yes
	$graph crosshairs configure -hide no
	$graph axis configure x -title "diameter (D), Angstroms"
	$graph axis configure y -title "$ext"
	pack append $graphWindow \
		$graph { fill expand padx 20 pady 20 }
	#
	# This window can be resized!
	wm min $graphWindow 10 10

	$graph postscript output $params(project)-$ext.ps \
	        -colormode color -width 5i -height 5i
        Blt_ZoomStack    $graph
        Blt_Crosshairs   $graph
        Blt_ActiveLegend $graph
        Blt_ClosestPoint $graph
	incr plotnum
}

######################################################################
#
#   read_input_file
#      reads the parameters from an input command file
#
proc read_input_file fileName {
  global params
  set f [open $fileName r]
  gets $f line; scan $line %s params(project)
  gets $f line; scan $line %s params(sasFile)
  gets $f line; scan $line %s%s params(qMin) params(qMax)
  gets $f line; scan $line %s params(contrast)
  gets $f line; scan $line %s params(dataFactor)
  gets $f line; scan $line %s params(esdFactor)
  gets $f line; scan $line %s params(bkg)
  gets $f line; scan $line %s params(shape)
  gets $f line; scan $line %s params(aspect)
  gets $f line; scan $line %s params(binType)
  gets $f line; scan $line %s params(nRadii)
  gets $f line; scan $line %s%s params(dMin) params(dMax)
  gets $f line; scan $line %s params(distType)
  gets $f line; scan $line %s params(defLevel)
  gets $f line; scan $line %s params(iterMax)
  gets $f line; scan $line %s params(slitLength)
  gets $f line; scan $line %s params(dE_E)
  gets $f line; scan $line %s params(method)
  close $f
}

######################################################################
#
#   write_input_file
#      writes the parameters to an input command file
#
proc write_input_file fileName {
  global params
  set f [open $fileName w]
  puts $f "[format %20s $params(project)] : Project Name  (only 1st item is read)"
  puts $f "[format %20s $params(sasFile)] : SAS file, contains columns: Q  i  esd"
  puts $f "[format %10s%10s $params(qMin) $params(qMax)] : qMin qMax, 1/A  (1.0e-8 to 100 means all data)"
  puts $f "[format %20s $params(contrast)] : rhosq       : scattering contrast, 10^20 1/cm^-4"
  puts $f "[format %20s $params(dataFactor)] : fac         :   I = fac * ( i - bkg )"
  puts $f "[format %20s $params(esdFactor)] : err         : ESD = fac * err * esd"
  puts $f "[format %20s $params(bkg)] : bkg         :   I = fac * ( i - bkg )"
  puts $f "[format %20s $params(shape)] : shapeModel  (1=spheroids, no others yet)"
  puts $f "[format %20s $params(aspect)] : Aspect Ratio"
  puts $f "[format %20s $params(binType)] : Bin Type    (1=Lin, 0=Log)"
  puts $f "[format %20s $params(nRadii)] : nRadii"
  puts $f "[format %10s%10s $params(dMin) $params(dMax)] : dMin dMax, A"
  puts $f "[format %20s $params(distType)] : n, in N(D)*V^n, 0=N(D), 1=f(D), 2=i(D)"
  puts $f "[format %20s $params(defLevel)] : defaultDistLevel  (MaxEnt only)"
  puts $f "[format %20s $params(iterMax)] : IterMax"
  puts $f "[format %20s $params(slitLength)] : slitLength, 1/A"
  puts $f "[format %20s $params(dE_E)] : dLambda/Lambda"
  puts $f "[format %20s $params(method)] : method (0=reg., 1=MaxEnt, 2=reg+NNLS, 3=NNLS, 4=SVD)"
  close $f
}

######################################################################
#
#   doOpen
#      reads the parameters to an input command file
#
proc doOpen {} {
	global params
	set oldDir [pwd]
	cd $params(wdName)
	read_input_file $params(cmdFile)
	cd $oldDir
}

######################################################################
#
#   doSave
#      saves the parameters to an input command file
#
proc doSave {} {
	global params
	if {![llength $params(project)]} return
	set oldDir [pwd]
	cd $params(wdName)
	write_input_file $params(cmdFile)
	cd $oldDir
}

######################################################################
#
#   doRun
#      runs the analysis
#
proc doRun {} {
	global params
	if {![llength $params(project)]} return
	set oldDir [pwd]
	cd $params(wdName)
	doSave
	_startAnalysis $params(project).log "sizes $params(cmdFile)"
	cd $oldDir
}

######################################################################
#
#   logWindow
#      writes into the log window
#
proc doViewFile s {
	global params
	if {$params(distType) == 0} {set ext "N-D"}
	if {$params(distType) == 1} {set ext "f-D"}
	if {$params(distType) == 2} {set ext "i-D"}
	if {$s == "cmd"} {set w .cmdw; set theFile $params(cmdFile)}
	if {$s == "SAS"} {set w .sasw; set theFile $params(sasFile)}
	if {$s == "log"} {set w .logw; set theFile $params(project).log}
	if {$s == "fit"} {set w .fitw; set theFile $params(project).fit}
	if {$s == "siz"} {set w .sizw; set theFile $params(project).$ext}
	if {$theFile == ""} return
	if {[winfo exists $w]} {destroy $w}
	toplevel $w
	wm title $w $theFile
	button $w.dismiss -text Dismiss 
	$w.dismiss configure -command [list destroy [winfo parent $w.dismiss]]
	pack $w.dismiss -side bottom -fill x
	text $w.text -yscrollcommand "$w.scroll set"
	scrollbar $w.scroll -command "$w.text yview"
	pack $w.scroll -side right -fill y
	pack $w.text -side left -fill both
	$w.text delete 1.0 end
	set f [open $theFile]
	while {![eof $f]} {
		$w.text insert end [read $f 1000]
	}
	close $f
	wm minsize $w 200 50
}

######################################################################
#
#   logWindow
#      writes into the log window
#
proc logWindow { text } {
    global params
    $params(textWidget) configure -state normal
    $params(textWidget) insert end $text\n
    $params(textWidget) see end
    $params(textWidget) configure -state disabled
}

######################################################################
#
#   Reader
#      reads command pipelines
#
proc Reader { pipe } {
  if [eof $pipe] {
    catch {close $pipe}
    logWindow "scan finished"
    changeButton scan
    catch {makeGnuplot}
    return
  }
  gets $pipe line
  logWindow $line
}

######################################################################
#
#   changeButton
#      for analyze or stop modes
#
proc changeButton { state } {
  global params pipe
  switch -- $state {
    analyze {
      bind . <Control-c> {}
      bind . <Control-a> doRun
    }
    stop {
      bind . <Control-c> "_stopAnalysis $pipe"
      bind . <Control-a> {}
    }
  }
}

######################################################################
#
#   _startAnalysis
#      starts an analysis running
#
proc _startAnalysis { logFile unixCommand } {
  global pipe
  logWindow "\n$logFile"
  logWindow "Started analysis [clock format [clock seconds]]"
  if [catch {
    set pipe [open "|$unixCommand |& cat" r]
    fileevent $pipe readable [list Reader $pipe]
  } result] {
    logWindow $result
  }
  changeButton stop
}

######################################################################
#
#   _stopAnalysis
#      stops a running analysis (assumes $pipe is an open pipe)
#
proc _stopAnalysis { pipe } {
  set pid [pid $pipe]
  exec kill [lindex [lsort $pid] 0]
  changeButton analysis
  logWindow "analysis interrupted"
}


######################################################################
#
#   grid_header
#      makes a 2-column header using the grid window manager
#
proc grid_header {widget text} {
  label $widget -text $text
  grid  configure $widget -columnspan 2 -sticky ew
}

######################################################################
#
#   grid_label_entry
#      makes a 2-column labeled entry using the grid window manager
#
proc grid_label_entry {widget label tclVar} {
  label [set widget]L -text $label
  entry [set widget]E -textvariable $tclVar -bg mintcream
  grid [set widget]L [set widget]E 
  grid configure [set widget]E -sticky ew
}

######################################################################
#
#   defineMenus
#      set the menus for the GUI
#
proc defineMenus {} {

  frame .mbar -relief raised -bd 2
  #frame .dummy -width 10c -height 1
  #pack .mbar .dummy -side top -fill x
  grid .mbar -columnspan 2 -sticky ew

  #
  #  menu File { New Open Save SaveAs Reload --- Quit }
  #
  menubutton .mbar.file -text File -menu .mbar.file.menu
  menu .mbar.file.menu
     .mbar.file.menu add command -label "Run the analysis" -command doRun
     .mbar.file.menu add separator
     .mbar.file.menu add command -label "Open command file" -command doOpen
     .mbar.file.menu add command -label "Save command file" -command doSave
     .mbar.file.menu add separator
     .mbar.file.menu add command -label "Quit" -command exit

  menubutton .mbar.plot -text Plot -menu .mbar.plot.menu
  menu .mbar.plot.menu
     .mbar.plot.menu add command -label SAS -command doSasPlot
     .mbar.plot.menu add command -label residuals -command doResPlot
     .mbar.plot.menu add command -label distribution -command doSizPlot

  menubutton .mbar.view -text Review -menu .mbar.view.menu
  menu .mbar.view.menu
     .mbar.view.menu add command -label "Command file" -command "doViewFile cmd"
     .mbar.view.menu add command -label "SAS file" -command "doViewFile SAS"
     .mbar.view.menu add separator
     .mbar.view.menu add command -label "Log file" -command "doViewFile log"
     .mbar.view.menu add command -label "Fit file" -command "doViewFile fit"
     .mbar.view.menu add command -label "Size file" -command "doViewFile siz"

  #grid 
  #grid configure .mbar.file -sticky w
  #grid configure .mbar.plot -sticky w
  #grid configure .mbar.view -sticky w
  pack .mbar.file .mbar.plot .mbar.view -side left

  tk_menuBar .mbar .mbar.file .mbar.edit .mbar.run  \
	.mbar.plot .mbar.view
}

######################################################################
#
#   define_GUI
#      define the GUI for this tool
#
proc define_GUI {} {
  global params

  grid_header .title "GUI for the \"sizes\" program" 
  frame .names   -relief raised -bd 2
  frame .basics  -relief raised -bd 2
  frame .misc    -relief raised -bd 2
  frame .options -relief raised -bd 2
  frame .text    -relief sunken -bd 2
  frame .buttons
    grid \
      [button .buttons.open  -text Open  -command doOpen] \
      [button .buttons.save  -text Save  -command doSave] \
      [button .buttons.run   -text Run   -command doRun] \
      [button .buttons.saxs  -text SAXS  -command doSasPlot] \
      [button .buttons.resid -text residuals -command doResPlot] \
      [button .buttons.dist  -text distribution  -command doSizPlot] \
      [button .buttons.quit  -text Quit  -command exit] \
    -sticky EW
  .buttons.quit configure -fg white -bg black \
    -activebackground yellow -activeforeground red

  #
  # work out a way to select either a configuration page 
  # or a text output page using a row of buttons at the top
  # of the window
  #
  # This will mean renaming the widgets and placing them into 
  # container frames
  #
  grid .names .options   -sticky nesw
  grid .basics .misc     -sticky nesw
  grid .text             -sticky ew  -columnspan 2
  grid .buttons          -sticky esw -columnspan 2

  grid_header .names.head "input file parameters"
  grid_label_entry .names.wd  "directory"     params(wdName)  
  grid_label_entry .names.cmd "command file"  params(cmdFile) 
  grid_label_entry .names.prj "project"       params(project) 
  grid_label_entry .names.sas "SAS data file" params(sasFile) 

  grid_header .basics.head "basic terms"
  grid_label_entry .basics.qMin   "Q_min, 1/A"    params(qMin) 
  grid_label_entry .basics.qMax   "Q_max, 1/A"    params(qMax) 
  grid_label_entry .basics.dMin   "D_min, 1/A"    params(dMin) 
  grid_label_entry .basics.dMax   "D_max, 1/A"    params(dMax) 
  grid_label_entry .basics.nRadii "\# of D bins"  params(nRadii) 
  grid_label_entry .basics.bkg    "background, input units" params(bkg) 
  label .basics.shapeL -text "Shape Model"
    set f .basics.shapeC; frame $f
    radiobutton $f.a -variable params(shape) -value 1 -text "spheroids"
    grid $f.a -sticky ew
    grid .basics.shapeL $f
    grid configure $f -sticky ew

  grid_header .misc.head "other factors"
  grid_label_entry .misc.fac      "data scale factor"    params(dataFactor) 
  grid_label_entry .misc.esd      "error scale factor"   params(esdFactor) 
  grid_label_entry .misc.contrast "contrast, 10^20/cm^4" params(contrast) 
  grid_label_entry .misc.aspect   "aspect ratio"         params(aspect) 
  grid_label_entry .misc.def      "default dist. level"  params(defLevel) 
  grid_label_entry .misc.dE_E     "dE_E"                 params(dE_E) 
  grid_label_entry .misc.slit     "slit length, 1/A"     params(slitLength) 

  grid_header .options.head options
  label .options.methodL -text "method"
    set f .options.methodC; frame $f
    grid \
      [radiobutton $f.me -variable params(method) -value 1 -text "MaxEnt"         ] \
      [radiobutton $f.re -variable params(method) -value 0 -text "Regularization" ] \
      [radiobutton $f.rn -variable params(method) -value 2 -text "Reg+NNLS"       ] \
      [radiobutton $f.nn -variable params(method) -value 3 -text "NNLS"           ] \
      [radiobutton $f.sv -variable params(method) -value 4 -text "SVD"            ] \
      -sticky ew
    grid .options.methodL $f
    grid configure $f -sticky ew
  label .options.binL -text "Binning interval"
    set f .options.binC; frame $f
    radiobutton $f.c -variable params(binType) -value 1 -text "Linear"
    radiobutton $f.d -variable params(binType) -value 0 -text "Log"
    grid $f.c $f.d -sticky ew
    grid .options.binL $f
    grid configure $f -sticky ew
  label .options.distL -text "Distribution"
    set f .options.distC; frame $f
    radiobutton $f.b -variable params(distType) -value 2 -text "i(D)"
    radiobutton $f.c -variable params(distType) -value 1 -text "f(D)"
    radiobutton $f.d -variable params(distType) -value 0 -text "N(D)"
    grid $f.b $f.c $f.d -sticky ew
    grid .options.distL $f
    grid configure $f -sticky ew
  scale .options.iterMax -label "Maximum iterations" \
                 -variable params(iterMax) \
                 -from 10 -to 300 -orient horizontal
  grid configure .options.iterMax -columnspan 2 -sticky ew

  # Create a text widget to log the output
  set f .text
    set    textCommand "text $f.log"
    append textCommand " -width 80 -height 15"
    append textCommand " -borderwidth 2"
    append textCommand " -relief ridge"
    append textCommand " -setgrid true"
    append textCommand " -yscrollcommand \{ $f.scroll set \}"
    append textCommand " -bg bisque"
    scrollbar $f.scroll -command [list $f.log yview]
    set params(textWidget) [eval $textCommand]
    grid $f.log $f.scroll -sticky nesw
    grid columnconfigure $f 0 -weight 1

  # the root is definitely not resizable
  # the text display window does not cooperate
  wm resizable . 0 0
}

######################################################################
#
#   setDefaults
#      set some defaults for each of the parameters
#
proc setDefaults {} {
  global params
  set defaults {
    qMin       1.0e-8   #  1/A
    qMax       100
    contrast   1        #  10^20/cm^4 (or 10^28/m^4)
    dataFactor 1        #    I = fac * ( i - bkg )
    esdFactor  1        #  ESD = fac * err * esd
    bkg        0.1      #    I = fac * ( i - bkg )
    shape      1        #  1=spheroids, no others yet
    aspect     1.0      #  r x r x r*aspect
    binType    1        #  1=Lin  0=Log
    nRadii     100
    dMin       25       #  Angstroms
    dMax       900
    distType   1        #  0=N(D)  1=f(D)  2=i(D)
    defLevel   1.0e-8
    iterMax    30       #  used only by MaxEnt
    slitLength 0        #  for slit-smeared data
    dE_E       0.0002   #  incident wavelength spread
    method     1        #  1=MaxEnt  0=regularization, 2=reg. with NNLS
    wdName     ./
    project    test
    sasFile    test.sas
    cmdFile    test.cmd
    font       *Times-Bold-R*14*   # for plot titles
  }
  foreach pair [split [string trim $defaults] \r\n] {
    set params([lindex $pair 0]) [lindex $pair 1]
  }
}

lappend auto_path $blt_library
set plotnum 1
setDefaults
defineMenus
define_GUI
changeButton analyze

if $argc {
  if {$argc > 1} {
    puts "usage:  tool [command_file]"
    exit
  }
  set params(cmdFile) $argv
  doOpen
}

