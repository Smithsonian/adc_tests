# This script freezes 'test_adc5g' and produces a binary

# Create the build directory if necessary
if [ ! -d "build" ]; then
    echo "A 'build' directory was not found, creating one."
    mkdir build
fi
cd build

# Find the Python freeze script, if it's not available exit
DEFAULTFREEZEPATH="/usr/share/doc/python2.6/examples/Tools/freeze/freeze.py"
if [ "$FREEZEPATH" == "" ]; then
    # FREEZEPATH not set, try default path
    echo "FREEZEPATH not set, looking for freeze.py in default location..."
    if [ ! -f $DEFAULTFREEZEPATH ]; then
	# Can't find freeze.py, GTFO
	echo "    Cannot find 'freeze.py', install python2.6-examples."
	echo "    If python2.6-examples is installed use 'dpkg -S freeze.py "
	echo "    and point FREEZEPATH to the location of 'freeze.py'"
	exit 1
    else
	echo "    Found 'freeze.py' in default location."
	FREEZELOC="$DEFAULTFREEZEPATH"
    fi
else
    FREEZELOC="$FREEZEPATH"
fi

# Use the freeze script to create a binary
python2.6 "$FREEZELOC" \
    -X _warnings -X copy \
    -X distutils -X locale -X macpath \
    -X ntpath -X os2emxpath -X popen2 -X pydoc \
    -X adc5g.mlab_tools \
    ../test_adc5g.py \
    -m encodings.ascii
make
cp test_adc5g ../bin/.