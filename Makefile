SAGE ?= sage

test:
	$(SAGE) -python -m pytest tests -q

doctest:
	PYTHONPATH=$(CURDIR) $(SAGE) -t diophantine_classifier/

coverage:
	$(SAGE) --coverage diophantine_classifier/

references:
	$(SAGE) -python tools/check_references.py

check: test doctest

smoke:
	$(SAGE) -python -c "from diophantine_classifier import classify; \
print(classify('x^2 - 61*y^2 = 1').explain())"

install:
	$(SAGE) -pip install -e .

.PHONY: test doctest coverage references check smoke install
