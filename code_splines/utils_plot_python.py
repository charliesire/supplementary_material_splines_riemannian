import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pyvista as pv
import pandas as pd
import seaborn as sns

def plot_sphere_3D(dfs, obs_points=None, obs_size=5, 
                                  sphere_opacity=1.0, n_rows=1, n_cols=None,
                                  titles=None, cmap=None,
                                  elev=30, azim=45, distance=10, save = False, file_name = ""):
    from matplotlib.colors import LinearSegmentedColormap
    if cmap is None:
        custom_colors = ["purple", "blue", "green", "white", "yellow", "orange", "red"]
        cmap = LinearSegmentedColormap.from_list("custom_cmap", custom_colors)

    dfs = [np.asarray(df) for df in dfs]
    n_spheres = len(dfs)

    if n_cols is None:
        n_cols = int(np.ceil(n_spheres / n_rows))

    all_values = np.concatenate([df[:, 2] for df in dfs])
    vmin, vmax = round(all_values.min(), 2), round(all_values.max(), 2)

    px_per_subplot = 1500  # haute qualité

    plotter = pv.Plotter(shape=(n_rows, n_cols), window_size=[px_per_subplot*n_cols, px_per_subplot*n_rows])
    plotter.set_background("white") 

    # Conversion elev/azim → position caméra 3D
    elev_rad = np.deg2rad(elev)
    azim_rad = np.deg2rad(azim)
    cam_x = distance * np.cos(elev_rad) * np.cos(azim_rad)
    cam_y = distance * np.cos(elev_rad) * np.sin(azim_rad)
    cam_z = distance * np.sin(elev_rad)
    camera_position = [(cam_x, cam_y, cam_z),  # position
                       (0, 0, 0),              # focal point
                       (0, 0, 1)]              # view-up

    for idx, df in enumerate(dfs):
        row = idx // n_cols
        col = idx % n_cols
        plotter.subplot(row, col)

        n_theta = len(np.unique(df[:, 0]))
        n_phi = len(np.unique(df[:, 1]))

        theta = np.concatenate([(df[:, 0] * np.pi / 180).reshape(n_phi, n_theta),
                                [df[df[:, 1] == 0, 0] * np.pi / 180]])
        phi = np.concatenate([(df[:, 1] * np.pi / 180).reshape(n_phi, n_theta),
                              [[0] * n_theta]])
        f = np.concatenate([(df[:, 2]).reshape(n_phi, n_theta),
                            [df[df[:, 1] == 0, 2]]])

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        grid = pv.StructuredGrid(x, y, z)
        grid['Test function'] = np.transpose(f).ravel()

        plotter.add_mesh(grid, scalars='Test function', cmap=cmap, 
                         clim=[vmin, vmax], opacity=sphere_opacity, label = None, scalar_bar_args= {"title" : "", "vertical": True, "position_x" : 0.86, "position_y" : 0.28, "bold":True, "label_font_size" : 50, "fmt": "%.2f","color": "black"})

        if obs_points is not None:
            obs_theta = np.deg2rad(obs_points[:, 0])
            obs_phi = np.deg2rad(obs_points[:, 1])
            ox = np.sin(obs_theta) * np.cos(obs_phi)
            oy = np.sin(obs_theta) * np.sin(obs_phi)
            oz = np.cos(obs_theta)

            points = np.column_stack((ox, oy, oz))
            plotter.add_points(points, color='black', point_size=obs_size,render_points_as_spheres=True)

        if titles and idx < len(titles):
            plotter.add_text(titles[idx], font_size=30, color="black")
        #plotter.show_axes()
        
        plotter.add_axes(line_width=4, xlabel='X', ylabel='Y', zlabel='Z', color='black',label_size=(1, 0.15))

        plotter.camera_position = camera_position
    if save: plotter.screenshot(file_name) 

    plotter.show()
    

def plot_sphere_2D(dfs, obs_points=None, n_rows=1, n_cols=None,
                            titles=None, cmap=None, save=False, file_name=""):
    """
    dfs: list of arrays [n_points, 3] with columns (theta[deg], phi[deg], value)
    theta in [0,180], phi in [0,360)
    """
    if cmap is None:
        custom_colors = ["purple", "blue", "green", "white", "yellow", "orange", "red"]
        cmap = LinearSegmentedColormap.from_list("custom_cmap", custom_colors)

    dfs = [np.asarray(df) for df in dfs]
    n_plots = len(dfs)

    if n_cols is None:
        n_cols = int(np.ceil(n_plots / n_rows))

    # Limites de couleur globales
    all_values = np.concatenate([df[:, 2] for df in dfs])
    vmin, vmax = np.nanmin(all_values), np.nanmax(all_values)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3.5*n_rows),
                             squeeze=False)

    for idx, df in enumerate(dfs):
        ax = axes[idx // n_cols, idx % n_cols]

        # Tri par theta puis phi
        df_sorted = df[np.lexsort((df[:, 1], df[:, 0]))]

        theta_unique = np.unique(df_sorted[:, 0])
        phi_unique = np.unique(df_sorted[:, 1])

        # Création de la grille
        TH, PH = np.meshgrid(theta_unique, phi_unique, indexing='ij')
        F = df_sorted[:, 2].reshape(len(theta_unique), len(phi_unique))

        pcm = ax.pcolormesh(PH, TH, F, cmap=cmap, shading='auto',
                            vmin=vmin, vmax=vmax)

        if obs_points is not None:
            ax.scatter(obs_points[:, 1], obs_points[:, 0],
                       c='black', s=20, marker='o',
                       edgecolors='white', linewidths=0.5)

        if titles and idx < len(titles):
            ax.set_title(titles[idx])

        ax.set_xlabel(r"$\phi$ [deg]")
        ax.set_ylabel(r"$\theta$ [deg]")
        ax.set_xlim(0, 360)
        ax.set_ylim(0, 180)
        ax.invert_yaxis()
    plt.tight_layout()
    fig.subplots_adjust(right=0.85)  

    cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
    fig.colorbar(pcm, cax=cbar_ax)




    if save:
        plt.savefig(file_name, dpi=600)
    plt.show()


def plot_cylinder_3D(dfs, obs_points=None, obs_size=5, 
                                    cylinder_opacity=1.0, n_rows=1, n_cols=None,
                                    titles=None, cmap=None, radius=1.0,
                                    elev=30, azim=45, distance=10,
                                    save=False, file_name=""):
    if cmap is None:
        custom_colors = ["purple", "blue", "green", "white", "yellow", "orange", "red"]
        cmap = LinearSegmentedColormap.from_list("custom_cmap", custom_colors)

    dfs = [np.asarray(df) for df in dfs]
    n_cylinders = len(dfs)

    if n_cols is None:
        n_cols = int(np.ceil(n_cylinders / n_rows))

    all_values = np.concatenate([df[:, 2] for df in dfs])
    vmin, vmax = all_values.min(), all_values.max()

    px_per_subplot = 1500  # haute qualité

    plotter = pv.Plotter(shape=(n_rows, n_cols), window_size=[px_per_subplot*n_cols, px_per_subplot*n_rows])
    plotter.set_background("white")

    elev_rad = np.deg2rad(elev)
    azim_rad = np.deg2rad(azim)
    cam_x = distance * np.cos(elev_rad) * np.cos(azim_rad)
    cam_y = distance * np.cos(elev_rad) * np.sin(azim_rad)
    cam_z = distance * np.sin(elev_rad)
    camera_position = [(cam_x, cam_y, cam_z), (0, 0, 0), (-1,0, 0)]


    for idx, df in enumerate(dfs):
        row = idx // n_cols
        col = idx % n_cols
        plotter.subplot(row, col)

        n_theta = len(np.unique(df[:, 0]))
        n_z = len(np.unique(df[:, 1]))

        # Convertir θ (en degrés) en radians
        theta = np.concatenate([(df[:, 0] * np.pi / 180).reshape(n_theta, n_z),
                                [df[df[:, 0] == 0, 0] * np.pi / 180]])
        z_vals = np.concatenate([df[:, 1].reshape(n_theta, n_z),
                                [df[df[:, 0] == 0, 1]]])
        f_vals = np.concatenate([df[:, 2].reshape(n_theta, n_z),
                                [df[df[:, 0] == 0, 2]]])

        # Conversion cylindre : rayon constant
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        z = z_vals - 1*np.mean(df[:,1])

        # Grille structurée pour PyVista
        grid = pv.StructuredGrid(x, y, z)
        grid['Test function'] = f_vals.ravel(order='F')

        plotter.add_mesh(grid, scalars='Test function', cmap=cmap,
                         clim=[vmin, vmax], opacity=cylinder_opacity,
                         scalar_bar_args={"title": "", "vertical": True,
                                          "position_x": 0.89, "position_y": 0.28,
                                          "bold": True, "label_font_size": 60, "fmt" :"%.1f","color": "black"})

        # Points d'observation
        if obs_points is not None:
            obs_theta = np.deg2rad(obs_points[:, 0])
            obs_z = obs_points[:, 1]
            ox = radius * np.cos(obs_theta)
            oy = radius * np.sin(obs_theta)
            oz = obs_z- 1*np.mean(df[:,1])
            points = np.column_stack((ox, oy, oz))
            plotter.add_points(points, color='black', point_size=obs_size,
                               render_points_as_spheres=True)

        if titles and idx < len(titles):
            plotter.add_text(titles[idx], font_size=40,color="black")

        plotter.add_axes(line_width=4, xlabel='X', ylabel='Y', zlabel='Z', color='black',label_size=(1, 0.15))

        
        plotter.camera_position = camera_position

    if save:
        plotter.screenshot(file_name)

    plotter.show()



def plot_cylinder_2D(dfs, obs_points=None, n_rows=1, n_cols=None,
                            titles=None, cmap=None, save=False, file_name=""):

    if cmap is None:
        custom_colors = ["purple", "blue", "green", "white", "yellow", "orange", "red"]
        cmap = LinearSegmentedColormap.from_list("custom_cmap", custom_colors)

    dfs = [np.asarray(df) for df in dfs]
    n_plots = len(dfs)

    if n_cols is None:
        n_cols = int(np.ceil(n_plots / n_rows))

    # Limites de couleur globales
    all_values = np.concatenate([df[:, 2] for df in dfs])
    vmin, vmax = np.nanmin(all_values), np.nanmax(all_values)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3.5*n_rows),
                             squeeze=False)

    for idx, df in enumerate(dfs):
        ax = axes[idx // n_cols, idx % n_cols]
        
        df_complete = df[df[:,0] == 0,]
        df_complete[:,0] = 360
        df = np.concatenate([df,df_complete], axis = 0)
        # Tri par theta puis phi
        df_sorted = df[np.lexsort((df[:, 1], df[:, 0]))]

        theta_unique = np.unique(df_sorted[:, 0])
        z_unique = np.unique(df_sorted[:, 1])

        # Création de la grille
        theta, z = np.meshgrid(theta_unique, z_unique, indexing='ij')
        F = df_sorted[:, 2].reshape(len(theta_unique), len(z_unique))

        pcm = ax.pcolormesh(z, theta, F, cmap=cmap, shading='auto',
                            vmin=vmin, vmax=vmax)

        if obs_points is not None:
            ax.scatter(obs_points[:, 1], obs_points[:, 0],
                       c='black', s=20, marker='o',
                       edgecolors='white', linewidths=0.5)

        if titles and idx < len(titles):
            ax.set_title(titles[idx])

        ax.set_xlabel(r"z")
        ax.set_ylabel(r"$\theta$ [deg]")
        ax.set_xlim(min(df[:,1]), max(df[:,1]))
        ax.set_ylim(min(df[:,0]), max(df[:,0]))
    plt.tight_layout()
    fig.subplots_adjust(right=0.88)  # Réduit la zone des subplots à gauche

    cbar_ax = fig.add_axes([0.9, 0.18, 0.03, 0.7])
    fig.colorbar(pcm, cax=cbar_ax)


    if save:
        plt.savefig(file_name, dpi=600)
    plt.show()

