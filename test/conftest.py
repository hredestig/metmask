import os
import pytest

from metmask.dbi import db, determine_path

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


@pytest.fixture(scope='session')
def memory_db():
    return db(':memory:')


@pytest.fixture(scope='session')
def package_db():
    db_file = os.path.join(determine_path(), 'data', 'metmask-db')
    return db(db_file)
