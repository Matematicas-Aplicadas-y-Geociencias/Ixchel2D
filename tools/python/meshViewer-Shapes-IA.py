import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.patches import Rectangle, Circle, Patch
import sys

class MeshShapeViewer:
    def __init__(self, x_coords, y_coords):
        self.x_coords = x_coords
        self.y_coords = y_coords
        self.shapes = []
        self.current_shape_type = None
        self.shape_params = {}
        self.properties_dict = {}  # Store properties for shapes
        self.property_counter = 0  # Counter for property IDs
        
        # Create figure and axes
        self.fig, self.ax = plt.subplots(figsize=(14, 10))
        plt.subplots_adjust(left=0.15, bottom=0.25)
        
        # Plot mesh
        self.plot_mesh()
        
        # Create UI elements
        self.create_buttons()
        self.create_textboxes()
        self.create_legend()
        
        # Bind mouse events
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        
    def plot_mesh(self):
        """Plot the mesh grid."""
        self.ax.clear()
        
        # Plot vertical lines
        for i in range(self.x_coords.shape[1]):
            self.ax.plot(self.x_coords[:, i], self.y_coords[:, i], 'b-', linewidth=0.5, alpha=0.6)
        
        # Plot horizontal lines
        for j in range(self.x_coords.shape[0]):
            self.ax.plot(self.x_coords[j, :], self.y_coords[j, :], 'b-', linewidth=0.5, alpha=0.6)
        
        # Plot points
        self.ax.scatter(self.x_coords, self.y_coords, c='blue', s=5, zorder=5, alpha=0.4)
        
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_title('2D Mesh with Interactive Shapes and Properties')
        self.ax.grid(True, alpha=0.2)
        self.ax.set_aspect('equal')
        
    def create_buttons(self):
        """Create control buttons."""
        button_props = {'color': 'lightblue', 'hovercolor': 'lightcyan'}
        
        # Rectangle button
        ax_rect_btn = plt.axes([0.12, 0.90, 0.08, 0.04])
        self.btn_rect = Button(ax_rect_btn, 'Rectangle', **button_props)
        self.btn_rect.on_clicked(self.select_rectangle)
        
        # Circle button
        ax_circle_btn = plt.axes([0.12, 0.85, 0.08, 0.04])
        self.btn_circle = Button(ax_circle_btn, 'Circle', **button_props)
        self.btn_circle.on_clicked(self.select_circle)
        
        # Clear last shape button
        ax_clear_btn = plt.axes([0.12, 0.80, 0.08, 0.04])
        self.btn_clear = Button(ax_clear_btn, 'Clear Last', color='lightyellow', hovercolor='lemonchiffon')
        self.btn_clear.on_clicked(self.clear_last_shape)
        
        # Clear all shapes button
        ax_clear_all_btn = plt.axes([0.12, 0.75, 0.08, 0.04])
        self.btn_clear_all = Button(ax_clear_all_btn, 'Clear All', color='lightcoral', hovercolor='salmon')
        self.btn_clear_all.on_clicked(self.clear_all_shapes)
        
    def create_textboxes(self):
        """Create input text boxes for shape parameters and properties."""
        # ===== SHAPE DEFINITION =====
        ax_rect_x = plt.axes([0.12, 0.68, 0.08, 0.03])
        self.textbox_rect_x = TextBox(ax_rect_x, 'Rect X:', initial='0.2', color='white')
        
        ax_rect_y = plt.axes([0.12, 0.64, 0.08, 0.03])
        self.textbox_rect_y = TextBox(ax_rect_y, 'Rect Y:', initial='0.2', color='white')
        
        ax_rect_w = plt.axes([0.12, 0.60, 0.08, 0.03])
        self.textbox_rect_w = TextBox(ax_rect_w, 'Rect W:', initial='0.3', color='white')
        
        ax_rect_h = plt.axes([0.12, 0.56, 0.08, 0.03])
        self.textbox_rect_h = TextBox(ax_rect_h, 'Rect H:', initial='0.2', color='white')
        
        # Draw rectangle button
        ax_draw_rect = plt.axes([0.12, 0.52, 0.08, 0.03])
        self.btn_draw_rect = Button(ax_draw_rect, 'Draw Rect', color='lightgreen', hovercolor='palegreen')
        self.btn_draw_rect.on_clicked(self.draw_rectangle)
        
        # Circle parameters
        ax_circle_x = plt.axes([0.12, 0.44, 0.08, 0.03])
        self.textbox_circle_x = TextBox(ax_circle_x, 'Circle X:', initial='0.5', color='white')
        
        ax_circle_y = plt.axes([0.12, 0.40, 0.08, 0.03])
        self.textbox_circle_y = TextBox(ax_circle_y, 'Circle Y:', initial='0.5', color='white')
        
        ax_circle_r = plt.axes([0.12, 0.36, 0.08, 0.03])
        self.textbox_circle_r = TextBox(ax_circle_r, 'Circle R:', initial='0.1', color='white')
        
        # Draw circle button
        ax_draw_circle = plt.axes([0.12, 0.32, 0.08, 0.03])
        self.btn_draw_circle = Button(ax_draw_circle, 'Draw Circle', color='lightgreen', hovercolor='palegreen')
        self.btn_draw_circle.on_clicked(self.draw_circle)
        
        # ===== PROPERTY DEFINITIONS =====
        # gamma_momen
        ax_prop1 = plt.axes([0.12, 0.24, 0.08, 0.03])
        self.textbox_gamma_momen = TextBox(ax_prop1, 'gamma_momen:', initial='0.0', color='white')
        
        # gamma_energ
        ax_prop2 = plt.axes([0.12, 0.20, 0.08, 0.03])
        self.textbox_gamma_energ = TextBox(ax_prop2, 'gamma_energ:', initial='0.0', color='white')
        
        # fuente_con_t
        ax_prop3 = plt.axes([0.12, 0.16, 0.08, 0.03])
        self.textbox_fuente_con_t = TextBox(ax_prop3, 'fuente_con_t:', initial='0.0', color='white')
        
        # fuente_lin_t
        ax_prop4 = plt.axes([0.12, 0.12, 0.08, 0.03])
        self.textbox_fuente_lin_t = TextBox(ax_prop4, 'fuente_lin_t:', initial='0.0', color='white')
        
        # Assign properties button
        ax_assign_prop = plt.axes([0.12, 0.08, 0.08, 0.03])
        self.btn_assign_prop = Button(ax_assign_prop, 'Assign Props', color='lightyellow', hovercolor='lemonchiffon')
        self.btn_assign_prop.on_clicked(self.assign_properties)
        
        # Export button
        ax_export = plt.axes([0.12, 0.04, 0.08, 0.03])
        self.btn_export = Button(ax_export, 'Export Data', color='lightcyan', hovercolor='lightblue')
        self.btn_export.on_clicked(self.export_to_file)
        
    def create_legend(self):
        """Create a legend showing instructions."""
        legend_text = (
            'WORKFLOW:\n'
            '1. Select shape (Rectangle/Circle)\n'
            '2. Enter shape parameters\n'
            '3. Click "Draw Rect" or "Draw Circle"\n'
            '4. Enter property values\n'
            '5. Click "Assign Props"\n'
            '6. Repeat steps 1-5 for more shapes\n'
            '7. Click "Export Data"\n\n'
            'Shape Parameters:\n'
            '  Rectangle: X, Y, Width, Height\n'
            '  Circle: X, Y, Radius\n\n'
            'Properties (all optional):\n'
            '  • gamma_momen\n'
            '  • gamma_energ\n'
            '  • fuente_con_t\n'
            '  • fuente_lin_t\n\n'
            'Output: i, j, property1, property2, ...'
        )
        self.fig.text(0.72, 0.12, legend_text, 
                     ha='left', fontsize=8,
                     bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8),
                     family='monospace')
        
    def select_rectangle(self, event):
        """Set current shape type to rectangle."""
        self.current_shape_type = 'rectangle'
        self.btn_rect.color = 'lightgreen'
        self.btn_circle.color = 'lightblue'
        self.fig.canvas.draw_idle()
        print("Rectangle mode selected. Enter parameters and click 'Draw Rect'")
        
    def select_circle(self, event):
        """Set current shape type to circle."""
        self.current_shape_type = 'circle'
        self.btn_circle.color = 'lightgreen'
        self.btn_rect.color = 'lightblue'
        self.fig.canvas.draw_idle()
        print("Circle mode selected. Enter parameters and click 'Draw Circle'")
        
    def draw_rectangle(self, event):
        """Draw a rectangle based on text box parameters."""
        try:
            x = float(self.textbox_rect_x.text)
            y = float(self.textbox_rect_y.text)
            width = float(self.textbox_rect_w.text)
            height = float(self.textbox_rect_h.text)
            
            rect = Rectangle((x, y), width, height, 
                           linewidth=2, edgecolor='red', 
                           facecolor='red', alpha=0.3, zorder=10)
            self.ax.add_patch(rect)
            
            shape_id = len(self.shapes)
            self.shapes.append({
                'type': 'rectangle',
                'patch': rect,
                'params': {'x': x, 'y': y, 'width': width, 'height': height},
                'id': shape_id,
                'properties': {
                    'gamma_momen': None,
                    'gamma_energ': None,
                    'fuente_con_t': None,
                    'fuente_lin_t': None
                }
            })
            
            print(f"Rectangle {shape_id} drawn: X={x}, Y={y}, W={width}, H={height}")
            self.fig.canvas.draw_idle()
            
        except ValueError:
            print("Error: Invalid rectangle parameters. Please enter valid numbers.")
            
    def draw_circle(self, event):
        """Draw a circle based on text box parameters."""
        try:
            x = float(self.textbox_circle_x.text)
            y = float(self.textbox_circle_y.text)
            radius = float(self.textbox_circle_r.text)
            
            circle = Circle((x, y), radius, 
                          linewidth=2, edgecolor='green', 
                          facecolor='green', alpha=0.3, zorder=10)
            self.ax.add_patch(circle)
            
            shape_id = len(self.shapes)
            self.shapes.append({
                'type': 'circle',
                'patch': circle,
                'params': {'x': x, 'y': y, 'radius': radius},
                'id': shape_id,
                'properties': {
                    'gamma_momen': None,
                    'gamma_energ': None,
                    'fuente_con_t': None,
                    'fuente_lin_t': None
                }
            })
            
            print(f"Circle {shape_id} drawn: X={x}, Y={y}, R={radius}")
            self.fig.canvas.draw_idle()
            
        except ValueError:
            print("Error: Invalid circle parameters. Please enter valid numbers.")
            
    def assign_properties(self, event):
        """Assign properties to the last drawn shape."""
        if not self.shapes:
            print("Error: No shapes to assign properties to. Draw a shape first.")
            return
        
        try:
            # Get values from textboxes
            gamma_momen = self.textbox_gamma_momen.text.strip()
            gamma_energ = self.textbox_gamma_energ.text.strip()
            fuente_con_t = self.textbox_fuente_con_t.text.strip()
            fuente_lin_t = self.textbox_fuente_lin_t.text.strip()
            
            # Validate all values are numbers
            properties_to_assign = {}
            
            if gamma_momen:
                try:
                    properties_to_assign['gamma_momen'] = float(gamma_momen)
                except ValueError:
                    print("Error: gamma_momen must be a number.")
                    return
            
            if gamma_energ:
                try:
                    properties_to_assign['gamma_energ'] = float(gamma_energ)
                except ValueError:
                    print("Error: gamma_energ must be a number.")
                    return
            
            if fuente_con_t:
                try:
                    properties_to_assign['fuente_con_t'] = float(fuente_con_t)
                except ValueError:
                    print("Error: fuente_con_t must be a number.")
                    return
            
            if fuente_lin_t:
                try:
                    properties_to_assign['fuente_lin_t'] = float(fuente_lin_t)
                except ValueError:
                    print("Error: fuente_lin_t must be a number.")
                    return
            
            if not properties_to_assign:
                print("Error: Enter at least one property value.")
                return
            
            # Assign to the last shape
            last_shape = self.shapes[-1]
            for prop_name, prop_value in properties_to_assign.items():
                last_shape['properties'][prop_name] = prop_value
            
            print(f"\nProperties assigned to {last_shape['type']} {last_shape['id']}:")
            for prop_name, prop_value in properties_to_assign.items():
                print(f"  {prop_name} = {prop_value}")
            print()
            
        except Exception as e:
            print(f"Error assigning properties: {e}")
            
    def clear_last_shape(self, event):
        """Remove the last drawn shape."""
        if self.shapes:
            shape_info = self.shapes.pop()
            shape_info['patch'].remove()
            print(f"Removed last {shape_info['type']}")
            self.fig.canvas.draw_idle()
        else:
            print("No shapes to clear.")
            
    def clear_all_shapes(self, event):
        """Remove all shapes."""
        for shape_info in self.shapes:
            shape_info['patch'].remove()
        self.shapes.clear()
        print("All shapes cleared.")
        self.fig.canvas.draw_idle()
        
    def on_click(self, event):
        """Handle mouse click events."""
        if event.inaxes != self.ax:
            return
        
        # Print coordinates where user clicked
        if event.button == 3:  # Right click
            print(f"Clicked at: X={event.xdata:.4f}, Y={event.ydata:.4f}")
    
    def point_in_rectangle(self, px, py, x, y, width, height):
        """Check if point (px, py) is inside rectangle."""
        return x <= px <= x + width and y <= py <= y + height
    
    def point_in_circle(self, px, py, cx, cy, radius):
        """Check if point (px, py) is inside circle."""
        return (px - cx)**2 + (py - cy)**2 <= radius**2
    
    def get_points_in_shapes(self):
        """Get all mesh points that fall inside defined shapes with their properties.
        Returns indices (i, j) where:
        - i runs horizontally (x direction, columns)
        - j runs vertically (y direction, rows)
        """
        points_with_properties = []
        
        # For each point in the mesh using indices
        # Now: j runs vertically (rows), i runs horizontally (columns)
        for j in range(self.x_coords.shape[0]):  # rows (y direction/vertical)
            for i in range(self.x_coords.shape[1]):  # columns (x direction/horizontal)
                px = self.x_coords[j, i]
                py = self.y_coords[j, i]
                
                # Check which shapes contain this point
                for shape in self.shapes:
                    is_inside = False
                    
                    if shape['type'] == 'rectangle':
                        params = shape['params']
                        is_inside = self.point_in_rectangle(
                            px, py, 
                            params['x'], params['y'], 
                            params['width'], params['height']
                        )
                    
                    elif shape['type'] == 'circle':
                        params = shape['params']
                        is_inside = self.point_in_circle(
                            px, py,
                            params['x'], params['y'],
                            params['radius']
                        )
                    
                    if is_inside:
                        # Collect all properties for this point from this shape
                        point_data = {'i': i, 'j': j}
                        for prop_name in ['gamma_momen', 'gamma_energ', 'fuente_con_t', 'fuente_lin_t']:
                            if shape['properties'][prop_name] is not None:
                                point_data[prop_name] = shape['properties'][prop_name]
                        
                        points_with_properties.append(point_data)
                        break  # Use the first shape that contains this point
        
        return points_with_properties
    
    def export_to_file(self, event):
        """Export mesh points inside shapes with properties to file."""
        points = self.get_points_in_shapes()
        
        if not points:
            print("Warning: No points found inside shapes with assigned properties.")
            return
        
        filename = 'mesh_properties.txt'
        
        try:
            with open(filename, 'w') as f:
                # Write header with total number of points
                f.write("# Mesh Points with Properties\n")
                f.write(f"# Total number of points: {len(points)}\n")
                f.write("# i (horizontal, columns), j (vertical, rows), gamma_momen, gamma_energ, fuente_con_t, fuente_lin_t\n")
                f.write("# " + "="*80 + "\n\n")
                
                # Write data
                for point in points:
                    i = point['i']+1
                    j = point['j']+1
                    gamma_momen = point.get('gamma_momen', '')
                    gamma_energ = point.get('gamma_energ', '')
                    fuente_con_t = point.get('fuente_con_t', '')
                    fuente_lin_t = point.get('fuente_lin_t', '')
                    
                    f.write(f"{i}  {j}  {gamma_momen}  {gamma_energ}  {fuente_con_t}  {fuente_lin_t}\n")
            
            print(f"\n{'='*70}")
            print(f"Successfully exported {len(points)} points to '{filename}'")
            print(f"{'='*70}\n")
            
            # Print summary
            print("Summary of exported points:")
            prop_counts = {'gamma_momen': 0, 'gamma_energ': 0, 'fuente_con_t': 0, 'fuente_lin_t': 0}
            for point in points:
                for prop in prop_counts:
                    if prop in point:
                        prop_counts[prop] += 1
            
            for prop, count in prop_counts.items():
                if count > 0:
                    print(f"  {prop}: {count} points")
            
        except Exception as e:
            print(f"Error exporting to file: {e}")
    
    def get_shapes_data(self):
        """Return list of all drawn shapes."""
        return self.shapes

def load_mesh(filename):
    """Load mesh coordinates from a text file."""
    try:
        data = np.loadtxt(filename)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception as e:
        print(f"Error loading mesh file: {e}")
        sys.exit(1)

def reshape_mesh(coordinates, nx, ny):
    """Reshape flat coordinate array into a 2D grid."""
    if len(coordinates) != nx * ny:
        raise ValueError(f"Expected {nx*ny} points but got {len(coordinates)}")
    
    x_coords = coordinates[:, 0].reshape(ny, nx)
    y_coords = coordinates[:, 1].reshape(ny, nx)
    
    return x_coords, y_coords

def main():
    """Main function to display interactive mesh and shape viewer."""
    
    if len(sys.argv) < 2:
        print("Usage: python mesh_viewer_shapes.py <mesh_file> [nx] [ny]")
        print("  mesh_file: path to the text file containing mesh coordinates")
        print("  nx: number of points in x direction (optional)")
        print("  ny: number of points in y direction (optional)")
        sys.exit(1)
    
    filename = sys.argv[1]
    coordinates = load_mesh(filename)
    
    # Auto-detect or use provided dimensions
    if len(sys.argv) >= 4:
        nx = int(sys.argv[2])
        ny = int(sys.argv[3])
    else:
        total_points = len(coordinates)
        nx = ny = int(np.sqrt(total_points))
        if nx * ny != total_points:
            print(f"Warning: {total_points} points detected.")
            print("Please provide nx and ny as arguments.")
            sys.exit(1)
    
    # Reshape coordinates
    x_coords, y_coords = reshape_mesh(coordinates, nx, ny)
    
    # Create viewer
    viewer = MeshShapeViewer(x_coords, y_coords)
    
    # Print welcome message
    print("\n" + "="*70)
    print("INTERACTIVE MESH SHAPE VIEWER WITH MULTIPLE PROPERTIES")
    print("="*70)
    print("\nIndex Convention:")
    print("  i: horizontal direction (x, columns) - ranges from 0 to nx-1")
    print("  j: vertical direction (y, rows) - ranges from 0 to ny-1")
    print("\nAvailable Properties:")
    print("  1. gamma_momen  - Momentum relaxation parameter")
    print("  2. gamma_energ  - Energy relaxation parameter")
    print("  3. fuente_con_t - Concentrated source term")
    print("  4. fuente_lin_t - Linear source term")
    print("\nWorkflow:")
    print("  1. Draw shapes (rectangles/circles)")
    print("  2. Enter property values (can leave blank if not needed)")
    print("  3. Click 'Assign Props' to assign to last shape")
    print("  4. Repeat for more shapes")
    print("  5. Click 'Export Data' to save mesh_properties.txt")
    print("\nOutput file format:")
    print("  i, j, gamma_momen, gamma_energ, fuente_con_t, fuente_lin_t")
    print("="*70 + "\n")
    
    plt.show()

if __name__ == '__main__':
    main()
