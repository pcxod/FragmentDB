CREATE TABLE Fragment (
	Id	INTEGER NOT NULL,
	tag	VARCHAR(255),
	Reference	TEXT,
	Name	TEXT,
	comment    TEXT,
	PRIMARY KEY(Id)
);


CREATE TABLE Atom (
	Id	INTEGER NOT NULL,
	FragmentId	INTEGER NOT NULL,
	tag	    VARCHAR(255),
	Name	VARCHAR(255),
	element	VARCHAR(2),
	x	FLOAT,
	y	FLOAT,
	z	FLOAT,
PRIMARY KEY(Id),
  FOREIGN KEY(FragmentId)
    REFERENCES Fragment(Id)
      ON DELETE CASCADE
      ON UPDATE NO ACTION);


CREATE TABLE Restraints (
  Id INTEGER  NOT NULL,
  FragmentId  INTEGER NOT NULL,
  ShelxName CHAR(4),
  Atoms TEXT,
PRIMARY KEY(Id),
  FOREIGN KEY(Id)
    REFERENCES Fragment(Id)
      ON DELETE CASCADE
      ON UPDATE NO ACTION);

CREATE INDEX Atom_FK ON Atom(FragmentId);
CREATE INDEX Restraint_FK ON Restraints(FragmentId);
CREATE INDEX Fragment_Name ON Atom(Name);
CREATE INDEX AtomId ON Atom(Id);




