#!/bin/bash
KSEQPP=kseq++
DEFAULT_PREFIX=/usr/local
TEMPLATE_PATH=./lib/pkgconfig/kseq++_pc.template
PKGCONFIG_FILE=./lib/pkgconfig/kseq++.pc
INSTALL="install -v"

PREFIX="$DEFAULT_PREFIX"
if [ -n "$1" ]; then
	PREFIX="$1"
else
	echo "No PREFIX provided. Assuming '$DEFAULT_PREFIX'"
fi

function generate_pcfile {
	prefix="$PREFIX"
	version=$(git describe)
	template=$(cat "$TEMPLATE_PATH")
	eval "echo \"$template\"" > "$PKGCONFIG_FILE"
}

function install_files {
	dest="$PREFIX/include/$KSEQPP"
	pc_dest="$PREFIX/lib/pkgconfig"
	src=./src
	mkdir -p "$dest"
	mkdir -p "$pc_dest"
	find "$src" -type f | xargs -I{} -- $INSTALL -v {} "$dest"
	$INSTALL "$PKGCONFIG_FILE" "$pc_dest"
}

generate_pcfile
install_files
