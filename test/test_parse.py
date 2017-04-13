from os import path
from metmask.parse import main as parse_main
from conftest import TEST_DATA_DIR


def test_csv(memory_db):
    importer = parse_main.importer(memory_db, 'simple',
                                   path.join(TEST_DATA_DIR, 'test-identifiers.txt'),
                                   'test')
    importer.parser.process()
    result = memory_db.simpleQuery('7732-18-5', 'cas', 'kegg', external=False)
    assert 'c00001' in result[0][0]
