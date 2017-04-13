from metmask.mask import mask

compounds = {'citrate': {'cas': '77-92-9',
                         'kegg': 'c00158'}}


def test_db_statistics(package_db, capsys):
    package_db.stats(more=True)
    out, err = capsys.readouterr()
    assert 'metmask db@' in out
    assert 'Known identifiers' in out
    assert 'Confidence codes' in out


def test_simple_query(package_db):
    result = package_db.simpleQuery(compounds['citrate']['cas'], 'cas', 'kegg')
    assert compounds['citrate']['kegg'] in result[0][0]
    result = package_db.simpleQuery(compounds['citrate']['cas'], 'cas', 'error')
    assert result[0][0] is None


def test_insert_delete_mask(memory_db):
    memory_db.createIdTable('table1')
    memory_db.createIdTable('table2')
    un = mask()
    un.append('table1', 'foo', 'strong', 'test', 'baz')
    un.append('table2', 'bar', 'strong', 'test', 'baz')
    memory_db.setMask(un)
    result = memory_db.simpleQuery('foo', 'table1', 'table2')
    assert 'bar' in result[0][0]
    memory_db.dropMask(un)
    result = memory_db.simpleQuery('foo', 'table1', 'table2')
    assert len(result[0][0]) == 0
