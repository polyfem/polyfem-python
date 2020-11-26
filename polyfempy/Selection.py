import json


class Selection:
    """Object used to select sideset and bodies"""

    def __init__(self):
        self.body_ids = []
        self.boundary_sidesets = []

    def select_body_with_sphere(self, id, center, radius):
        """Select a body using a sphere"""
        self.body_ids.append({"id": id, "center": center, "radius": radius})

    def select_body_with_box(self, id, box_min, box_max):
        """Select a body using an axis-aligned box"""
        self.body_ids.append({"id": id, "box": [box_min, box_max]})

    def select_body_with_axis_plane(self, id, axis, position):
        """Select a body using an axis-aligned plane at position position, axis 1, 2, 3 is x, y, and z respectively. Use negative to flip axis (e.g., -1 is negative x axis)"""
        self.body_ids.append({"id": id, "position": position, "axis": axis})

    def select_body_with_plane(self, id, normal, offset):
        """Select a body using a generic plane with normal normal, the point on the plane is defined by normal*offset"""
        self.body_ids.append({"id": id, "normal": normal, "offset": offset})

    def select_sideset_with_sphere(self, id, center, radius):
        """Select a boundary sideset using a sphere"""
        self.boundary_sidesets.append(
            {"id": id, "center": center, "radius": radius})

    def select_sideset_with_box(self, id, box_min, box_max):
        """Select a boundary sideset using an axis-aligned box"""
        self.boundary_sidesets.append({"id": id, "box": [box_min, box_max]})

    def select_sideset_with_axis_plane(self, id, axis, position):
        """Select a boundary sideset using an axis-aligned plane at position position, axis 1, 2, 3 is x, y, and z respectively. Use negative to flip axis (e.g., -1 is negative x axis)"""
        self.boundary_sidesets.append(
            {"id": id, "position": position, "axis": axis})

    def select_sideset_with_plane(self, id, normal, offset):
        """Select a boundary sideset using a generic plane with normal normal, the point on the plane is defined by normal*offset"""
        self.boundary_sidesets.append(
            {"id": id, "normal": normal, "offset": offset})

    def __str__(self):
        """stringygied json description of this class, used to run the solver"""

        tmp = dict(
            (key, value)
            for (key, value) in self.__dict__.items())

        return json.dumps(tmp, sort_keys=True, indent=4)
