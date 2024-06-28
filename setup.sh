#!/usr/bin/env bash

OpenMMScripts_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

add_to_bashrc() {
    echo "" >> ~/.bashrc
    echo "# Automatically added by setup.sh script for OpenMMScripts repository" >> ~/.bashrc
    echo "export PATH=\"$OpenMMScripts_path:\$PATH\"" >> ~/.bashrc
}

confirm() {
    while true; do
        read -p "Do you want to add \"$OpenMMScripts_path\" to your PATH in .bashrc? (y/n): " yn
        case $yn in
            [Yy]* ) add_to_bashrc; break;;
            [Nn]* ) echo "Skipping modification of .bashrc."; break;;
            * ) echo "Please answer yes or no.";;
        esac
    done
}

confirm
