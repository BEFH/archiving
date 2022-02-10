import sqlalchemy as sa

engine = sa.create_engine('sqlite:////sc/arion/projects/LOAD/archive/archive.sqlite')

metadata_obj = sa.MetaData()

archive = sa.Table('archive', metadata_obj,
    sa.Column('archive_id', sa.String(16), primary_key=True),
    sa.Column('user_name', sa.String(16), nullable=False),
    sa.Column('time', sa.DateTime(timezone=True)),
    sa.Column('archive_name', sa.String(50), nullable=False),
    sa.Column('archive_directory', sa.String(192), nullable=False),
    sa.Column('total_mib', sa.Float),
    sa.Column('freed_mib', sa.Float),
    sa.Column('kept_mib', sa.Float)
)

file = sa.Table('file', metadata_obj,
    sa.Column('file_id', sa.String(16), primary_key=True),
    sa.Column('archive_id', sa.String(16), sa.ForeignKey("archive.archive_id"), nullable=False),
    sa.Column('filename', sa.String(64), nullable=False),
    sa.Column('path', sa.String(256), nullable=False),
    sa.Column('extension', sa.String(15)),
    sa.Column('kind', sa.String(15)),
    sa.Column('keep', sa.String(15)),
    sa.Column('directory', sa.Boolean),
    sa.Column('kind', sa.String(15)),
    sa.Column('size_mib', sa.Float),
    sa.Column('mode', sa.Float),
    sa.Column('user', sa.String(16)),
    sa.Column('group', sa.String(16)),
    sa.Column('time_modified', sa.DateTime(timezone=True)),
    sa.Column('removal', sa.String(16)),
)

metadata_obj.create_all(engine)

