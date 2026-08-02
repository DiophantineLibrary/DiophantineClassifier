SAGE ?= sage

test:
	$(SAGE) -python -m pytest tests -q

smoke:
	$(SAGE) -python -c "from diophantine_classifier import classify; \
print(classify('x^2 - 61*y^2 = 1').explain())"

install:
	$(SAGE) -pip install -e .

.PHONY: test smoke install
