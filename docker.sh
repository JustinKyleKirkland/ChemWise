#!/bin/bash

# ChemWise Docker Management Script
# This script provides convenient commands for building and running ChemWise in Docker

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Functions
print_help() {
    echo -e "${BLUE}ChemWise Docker Management Script${NC}"
    echo ""
    echo "Usage: $0 [COMMAND] [OPTIONS]"
    echo ""
    echo "Commands:"
    echo "  build              Build the ChemWise Docker image"
    echo "  build-jupyter      Build the ChemWise Jupyter Docker image"
    echo "  run                Run ChemWise interactively"
    echo "  run-example        Run an example calculation"
    echo "  jupyter            Start Jupyter Lab server"
    echo "  shell              Start an interactive shell in the container"
    echo "  clean              Remove all ChemWise Docker images and containers"
    echo "  logs               Show container logs"
    echo "  stop               Stop all running ChemWise containers"
    echo ""
    echo "Examples:"
    echo "  $0 build                    # Build the main image"
    echo "  $0 run water.xyz            # Run calculation on water.xyz"
    echo "  $0 run-example water        # Run the water example"
    echo "  $0 jupyter                  # Start Jupyter Lab"
    echo "  $0 shell                    # Interactive shell"
    echo ""
}

build_image() {
    echo -e "${GREEN}Building ChemWise Docker image...${NC}"
    docker build -t chemwise:latest .
    echo -e "${GREEN}Build completed successfully!${NC}"
}

build_jupyter() {
    echo -e "${GREEN}Building ChemWise Jupyter Docker image...${NC}"
    # First ensure main image exists
    if ! docker image inspect chemwise:latest >/dev/null 2>&1; then
        echo -e "${YELLOW}Main ChemWise image not found. Building it first...${NC}"
        build_image
    fi
    docker build -f Dockerfile.jupyter -t chemwise-jupyter:latest .
    echo -e "${GREEN}Jupyter build completed successfully!${NC}"
}

run_interactive() {
    echo -e "${GREEN}Starting ChemWise interactive session...${NC}"
    docker run -it --rm \
        -v "$SCRIPT_DIR/data:/app/data" \
        -v "$SCRIPT_DIR/output:/app/output" \
        -v "$SCRIPT_DIR/examples:/app/examples" \
        chemwise:latest \
        /bin/bash
}

run_calculation() {
    local input_file="$1"
    if [ -z "$input_file" ]; then
        echo -e "${RED}Error: Please provide an input file${NC}"
        echo "Usage: $0 run <input_file>"
        exit 1
    fi
    
    echo -e "${GREEN}Running ChemWise calculation on $input_file...${NC}"
    docker run --rm \
        -v "$SCRIPT_DIR/data:/app/data" \
        -v "$SCRIPT_DIR/output:/app/output" \
        -v "$SCRIPT_DIR/examples:/app/examples" \
        chemwise:latest \
        chemwise calculate "/app/data/$input_file"
}

run_example() {
    local example="$1"
    if [ -z "$example" ]; then
        echo -e "${RED}Error: Please provide an example name${NC}"
        echo "Available examples: water, methane, h2"
        exit 1
    fi
    
    echo -e "${GREEN}Running ChemWise example: $example${NC}"
    docker run --rm \
        -v "$SCRIPT_DIR/output:/app/output" \
        chemwise:latest \
        chemwise calculate "/app/examples/${example}.xyz"
}

start_jupyter() {
    echo -e "${GREEN}Starting Jupyter Lab server...${NC}"
    
    # Create notebooks directory if it doesn't exist
    mkdir -p "$SCRIPT_DIR/notebooks"
    
    docker run -d \
        --name chemwise-jupyter \
        -p 8888:8888 \
        -v "$SCRIPT_DIR/data:/app/data" \
        -v "$SCRIPT_DIR/output:/app/output" \
        -v "$SCRIPT_DIR/examples:/app/examples" \
        -v "$SCRIPT_DIR/notebooks:/app/notebooks" \
        chemwise-jupyter:latest
    
    echo -e "${GREEN}Jupyter Lab started!${NC}"
    echo -e "${BLUE}Access it at: http://localhost:8888${NC}"
    echo -e "${BLUE}Token: chemwise-token${NC}"
}

show_logs() {
    local service="$1"
    if [ -z "$service" ]; then
        service="chemwise-jupyter"
    fi
    
    echo -e "${GREEN}Showing logs for $service...${NC}"
    docker logs -f "$service"
}

stop_containers() {
    echo -e "${GREEN}Stopping ChemWise containers...${NC}"
    docker stop chemwise-jupyter 2>/dev/null || true
    docker stop chemwise-app 2>/dev/null || true
    echo -e "${GREEN}Containers stopped.${NC}"
}

clean_all() {
    echo -e "${YELLOW}Cleaning up ChemWise Docker resources...${NC}"
    
    # Stop containers
    stop_containers
    
    # Remove containers
    docker rm chemwise-jupyter 2>/dev/null || true
    docker rm chemwise-app 2>/dev/null || true
    
    # Remove images
    docker rmi chemwise:latest 2>/dev/null || true
    docker rmi chemwise-jupyter:latest 2>/dev/null || true
    
    echo -e "${GREEN}Cleanup completed.${NC}"
}

# Main script logic
case "$1" in
    "build")
        build_image
        ;;
    "build-jupyter")
        build_jupyter
        ;;
    "run")
        if [ "$2" ]; then
            run_calculation "$2"
        else
            run_interactive
        fi
        ;;
    "run-example")
        run_example "$2"
        ;;
    "jupyter")
        if ! docker image inspect chemwise-jupyter:latest >/dev/null 2>&1; then
            echo -e "${YELLOW}Jupyter image not found. Building it first...${NC}"
            build_jupyter
        fi
        start_jupyter
        ;;
    "shell")
        run_interactive
        ;;
    "logs")
        show_logs "$2"
        ;;
    "stop")
        stop_containers
        ;;
    "clean")
        clean_all
        ;;
    "help"|"--help"|"-h"|"")
        print_help
        ;;
    *)
        echo -e "${RED}Unknown command: $1${NC}"
        print_help
        exit 1
        ;;
esac
