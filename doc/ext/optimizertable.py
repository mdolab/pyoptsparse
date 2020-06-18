from importlib import import_module
from docutils.parsers.rst import directives
from docutils.parsers.rst.directives.tables import Table
from docutils.statemachine import ViewList
from docutils import nodes


class OptimizerTable(Table):
    """
    This Table directive formats the defOpts dictionary in a nice table

    This is heavily adapted from
    https://github.com/BU-NU-CLOUD-SP16/Trusted-Platform-Module-nova/blob/master/api-ref/ext/rest_parameters.py
    """

    required_arguments = 1
    has_content = True
    option_spec = {
        "type": directives.uri,
        "header": str,
        "widths": directives.positive_int_list,
    }

    def get_options_informs(self, cls):
        optimizer_instance = cls()
        self.defOpts = optimizer_instance.defOpts
        self.informs = optimizer_instance.informs

    def set_header(self):
        if "header" in self.options:
            self.header = self.options["header"].split(",")
        else:
            if self.options["type"] == "options":
                self.header = ["Option name", "Type", "Default value"]
            elif self.options["type"] == "informs":
                self.header = ["Code", "Description"]

        self.max_cols = len(self.header)

    def set_width(self):
        if "widths" in self.options:
            self.col_widths = self.options["widths"]
        else:
            if self.options["type"] == "options":
                self.col_widths = [20, 10, 30]
            elif self.options["type"] == "informs":
                self.col_widths = [5, 55]

    def run(self):
        module_path, member_name = self.arguments[0].rsplit(".", 1)
        class_name = getattr(import_module(module_path), member_name)
        # set the self.defOpts and self.informs attribute
        self.get_options_informs(class_name)
        # set header option
        self.set_header()
        # set width
        self.set_width()
        table_node = self.build_table()
        return [table_node]

    def collect_rows(self):
        # Add a column for a field. In order to have the RST inside
        # these fields get rendered, we need to use the
        # ViewList. Note, ViewList expects a list of lines, so chunk
        # up our content as a list to make it happy.
        def add_col(value):
            entry = nodes.entry()
            result = ViewList(value.split("\n"))
            self.state.nested_parse(result, 0, entry)
            return entry

        rows = []
        groups = []
        # options
        if self.options["type"] == "options":
            for key, values in self.defOpts.items():
                trow = nodes.row()
                # first add the name column, with text = key
                trow += add_col("``" + key + "``")
                # loop over dictionary values
                for idx, value in enumerate(values):
                    # this is the type of the default value, so we need to extract the __name__ attribute
                    if idx == 0:
                        trow += add_col(value.__name__)
                    # this is the default value
                    else:
                        if isinstance(value, (float, int)):
                            limit = 5
                            if abs(value) > 10 ** limit or (abs(value) < 10 ** (-limit) and abs(value) > 0):
                                value = "{:.3E}".format(value)
                            else:
                                value = "{}".format(value)
                        trow += add_col(str(value))
                rows.append(trow)
        # informs
        elif self.options["type"] == "informs":
            for key, value in self.informs.items():
                trow = nodes.row()
                # first add the name column, with text = key
                trow += add_col("``" + str(key) + "``")
                # add inform description
                trow += add_col(value)
                rows.append(trow)
        return rows, groups

    def build_table(self):
        table = nodes.table()
        tgroup = nodes.tgroup(cols=len(self.header))
        table += tgroup

        tgroup.extend(
            nodes.colspec(colwidth=col_width, colname="c" + str(idx)) for idx, col_width in enumerate(self.col_widths)
        )

        thead = nodes.thead()
        tgroup += thead

        row_node = nodes.row()
        thead += row_node
        row_node.extend(nodes.entry(h, nodes.paragraph(text=h)) for h in self.header)

        tbody = nodes.tbody()
        tgroup += tbody

        rows, groups = self.collect_rows()
        tbody.extend(rows)
        table.extend(groups)

        return table


def setup(app):
    app.add_directive("optimizertable", OptimizerTable)
