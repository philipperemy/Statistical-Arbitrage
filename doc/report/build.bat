@echo off
echo RUNNING...
:loop
	pdflatex Thesis.tex -synctex=1 -interaction=nonstopmode -shell-escape -output-directory out > nul
	move out\Thesis.pdf . > nul
	rem del /Q out\*
	ping -n 5 127.0.0.1 > nul
goto loop

