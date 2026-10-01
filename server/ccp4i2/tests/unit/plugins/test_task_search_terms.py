from ccp4i2.core.tasks import TASKS


def test_search_terms_are_tuples_of_lowercase_strings():
    for name, task in TASKS.items():
        assert isinstance(task.searchTerms, tuple), name
        for term in task.searchTerms:
            assert isinstance(term, str) and term == term.lower().strip(), (name, term)


def test_search_terms_add_words_the_task_lacks():
    for name, task in TASKS.items():
        fields = (name, task.title, task.shortTitle, task.description)
        text = " ".join(f or "" for f in fields).lower()
        for term in task.searchTerms:
            assert term not in text, (name, term)


def test_no_search_term_is_inside_another():
    for name, task in TASKS.items():
        for term in task.searchTerms:
            others = [t for t in task.searchTerms if t != term]
            assert not any(term in t for t in others), (name, term)
