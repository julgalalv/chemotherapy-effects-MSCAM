# chemotherapy-effects-SCAM
A cellular automata model of chemotherapy effects on tumour growth (based on https://doi.org/10.1080/13873954.2019.1571515)

```flow
st=>start: Start
op=>operation: Your Operation
cond=>condition: Yes or No?
e=>end

st->op->cond
cond(yes)->e
cond(no)->op
​```