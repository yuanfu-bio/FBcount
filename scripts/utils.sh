#! /bin/bash
# Generic Functions

# define some logger functions
timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

log_info() {
    echo "$(timestamp) [INFO] $1"
}

log_warning() {
    echo "$(timestamp) [WARNING] $1"
}

log_error() {
    echo "$(timestamp) [ERROR] $1"
}
