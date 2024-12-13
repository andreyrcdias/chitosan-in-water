run:
	./run.sh

run-v2:
	./run_v2.sh

run-v3:
	./run_v3.sh

clean:
	./clean.sh

init:
	sudo unzip -o "LCQC_oplsaa.ff.zip" -d "/usr/local/gromacs/share/gromacs/top/oplsaa.ff"
