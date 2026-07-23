SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-1 -t ghcr.io/pereiramemo/mg-clust/module-1:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-2 -t ghcr.io/pereiramemo/mg-clust/module-2:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-3 -t ghcr.io/pereiramemo/mg-clust/module-3:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-4 -t ghcr.io/pereiramemo/mg-clust/module-4:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-5 -t ghcr.io/pereiramemo/mg-clust/module-5:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-6 -t ghcr.io/pereiramemo/mg-clust/module-6:latest .