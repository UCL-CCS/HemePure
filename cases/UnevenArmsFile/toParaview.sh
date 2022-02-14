#!/usr/bin/env bash

processing () {

	NAME=${1?No file provided}
	OUT=${2?No output name provided}

	echo Processing $NAME for ParaView, to be called $OUT.

	sed -i '1,2d' $1
	sed -i '/^$/d' $1

	#awk -F ' ' 'BEGIN{N=-1; last_step=-1}{if(last_step != $1) {N++; print>$OUT"N".txt"; last_step=$1;} else {print>>$OUT"N".txt"}}' $NAME
	awk -v out=$2 -F ' ' 'BEGIN{N=-2; last_step=-1}{if(last_step != $1) {N++; print>out N".txt"; last_step=$1;} else {print>>out N".txt"}}' $NAME

	sed -i '1isteps gridX gridY gridZ velX velY velZ pressure' $OUT*
}

processing "results/outlet.txt" "results/outletPV"
processing "results/inlet.txt"  "results/inletPV"
processing "results/planeY.txt"  "results/planeYPV"

