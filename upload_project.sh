if [ -d dist ]
then
echo "Removing old dist directory"
rm -rf dist/
fi

echo "Building..."
python3 -m build

echo "Uploading... requires pypi user name and pw."
python3 -m twine upload dist/*

