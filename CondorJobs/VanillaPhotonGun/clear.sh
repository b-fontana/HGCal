if [ "${1}" == "yes" ]; then
    #rm output/*;
    #rm log/*;
    #rm error/*;
    echo "Files deleted.";
else
    echo "No files deleted. Use 'bash clear.sh yes' to delete them."
fi
