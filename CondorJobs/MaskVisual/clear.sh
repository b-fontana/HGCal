if [ "${1}" == "yes" ]; then
    rm out/*;
    rm log/*;
    echo "Files deleted.";
else
    echo "No files deleted. Use 'bash clear.sh yes' to delete them."
fi
