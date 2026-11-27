PYTHON ?= python
LATEX ?= pdflatex
.PHONY: check pdf clean
check:
	$(PYTHON) scripts/audit.py
	$(PYTHON) scripts/audit_good_turing.py
	$(PYTHON) scripts/audit_good_turing_phase.py
pdf:
	cd paper && $(LATEX) -interaction=nonstopmode -halt-on-error main.tex > build-pass1.log
	cd paper && $(LATEX) -interaction=nonstopmode -halt-on-error main.tex > build-pass2.log
	cd paper && $(LATEX) -interaction=nonstopmode -halt-on-error main.tex > build-pass3.log
	cd good-turing && $(LATEX) -interaction=nonstopmode -halt-on-error main.tex > build-pass1.log
	cd good-turing && $(LATEX) -interaction=nonstopmode -halt-on-error main.tex > build-pass2.log
	cd good-turing && $(LATEX) -interaction=nonstopmode -halt-on-error main.tex > build-pass3.log
clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc good-turing/*.aux good-turing/*.log good-turing/*.out good-turing/*.toc
