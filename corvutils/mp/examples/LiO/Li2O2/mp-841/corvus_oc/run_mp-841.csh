#!/bin/csh

# Remove the pending flag
rm CTRL_PENDING

# Set a flag so crv_mp_stat knows this is running
touch CTRL_RUNNING


# Set the Corvus location
setenv CORVUS_EXE /home/fdv/.local/bin/run-corvus

# Run Corvus
$CORVUS_EXE -i opcons.in

# Remove the running flag
rm CTRL_RUNNING

# Set a flag indicating the calculation is done
touch CTRL_COMPLETED

