makeEditDist:
	tar -xvzf packages/py-editdist-0.3.tar.gz
	cd py-editdist-0.3; python setup.py build
	mkdir -p resources/usr
	cd py-editdist-0.3; python setup.py install --prefix=../resources/usr

clean:
	rm -rf py-editdist-0.3
	rm -rf resources/usr
