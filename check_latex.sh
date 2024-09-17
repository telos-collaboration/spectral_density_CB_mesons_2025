#!/bin/bash

# Function to convert line endings
convert_line_endings() {
    local file=$1
    dos2unix "$file"
}

# Check if LaTeX is installed
if ! command -v latex > /dev/null 2>&1; then
    echo "LaTeX is not installed."

    # File to be modified
    FILE="./plateaus/Lib/plot_package.py"

    # Check if the file exists
    if [ -f "$FILE" ]; then
        # Backup the original file
        cp "$FILE" "$FILE.bak"
        echo "Backup of $FILE created as $FILE.bak"

        # Add the line to the end of the file
        echo "rcParams['text.usetex'] = False" >> "$FILE"
        echo "Modified $FILE to include rcParams['text.usetex'] = False"
    else
        echo "$FILE does not exist."
    fi

    # Find and modify all paperdraft.mplstyle files
    find . -type f -name 'paperdraft.mplstyle' | while read -r stylefile; do
        # Backup the original style file
        cp "$stylefile" "$stylefile.bak"
        echo "Backup of $stylefile created as $stylefile.bak"

        # Convert line endings
        convert_line_endings "$stylefile"
        
        # Add the line to the style file if it doesn't already contain it
        if ! grep -q 'text.usetex: False' "$stylefile"; then
            echo "text.usetex: False" >> "$stylefile"
            echo "Modified $stylefile to include text.usetex: False"
        else
            echo "$stylefile already contains text.usetex: False"
        fi
    done
else
    echo "LaTeX is installed."
fi
