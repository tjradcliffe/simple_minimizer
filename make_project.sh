if [ -d dist ]
then
echo "Removing old dist directory"
rm -rf dist/
fi

echo "Building..."
python3 -m build

echo "Uploading..."
python3 -m twine upload --repository testpypi dist/*

echo "Installing..."
python3 -m pip install --upgrade --index-url https://test.pypi.org/simple/ --no-deps simple-minimizer
