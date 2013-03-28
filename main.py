import sys

from search_gui import Dialog
from PySide.QtGui import QApplication

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = Dialog()
    dialog.show()
    dialog.initSearchDb()
    sys.exit(app.exec_())
