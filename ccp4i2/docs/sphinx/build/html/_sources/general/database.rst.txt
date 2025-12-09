########################
CCP4i2 — CCP4i2 Database
########################

::

    cd $HOME/.CCP4I2/db
    echo '.dump' | sqlite3 database.sqlite  > database.dump
    echo '.read database.dump' | sqlite rebuilt_database.sqlite

.. |CCP4 web logo| image:: ccp4i2.png
   :width: 75px
