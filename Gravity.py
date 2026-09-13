import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
import os

G = 6.67430e-11
AU = 1.495978707e11
DAY = 86400
MSUN = 1.98847e30

os.makedirs("outputs", exist_ok=True)


def get_number(prompt, minimum=None):
    while True:
        try:
            value = float(input(prompt))

            if minimum is not None and value < minimum:
                print(f"Enter a value >= {minimum}.")
                continue

            return value

        except ValueError:
            print("Enter a valid number.")


def position(radius_au, angle_deg):
    angle = np.radians(angle_deg)

    return np.array([
        radius_au * AU * np.cos(angle),
        radius_au * AU * np.sin(angle),
        0.0
    ])


def velocity(speed_kms, direction_deg):
    angle = np.radians(direction_deg)

    return np.array([
        speed_kms * 1000 * np.cos(angle),
        speed_kms * 1000 * np.sin(angle),
        0.0
    ])


def acceleration(positions, masses, softening=0.0):
    n = len(masses)

    acc = np.zeros_like(positions)

    for i in range(n):
        displacement = positions - positions[i]

        distance_squared = np.sum(
            displacement**2,
            axis=1
        )

        if softening > 0:
            distance_squared += softening**2

        distance_squared[i] = np.inf

        inv_distance_cubed = 1.0 / (
            distance_squared *
            np.sqrt(distance_squared)
        )

        acc[i] = G * np.sum(
            masses[:, None]
            * displacement
            * inv_distance_cubed[:, None],
            axis=0
        )

    return acc


def verlet_step(
    positions,
    velocities,
    masses,
    dt,
    softening=0.0
):
    acc = acceleration(
        positions,
        masses,
        softening
    )

    new_positions = (
        positions
        + velocities * dt
        + 0.5 * acc * dt**2
    )

    new_acc = acceleration(
        new_positions,
        masses,
        softening
    )

    new_velocities = (
        velocities
        + 0.5 * (acc + new_acc) * dt
    )

    return new_positions, new_velocities


def simulate(
    positions,
    velocities,
    masses,
    duration_days,
    dt_days,
    softening=0.0
):
    steps = int(duration_days / dt_days) + 1

    history = np.zeros(
        (steps, len(masses), 3)
    )

    history[0] = positions

    current_positions = positions.copy()
    current_velocities = velocities.copy()

    dt = dt_days * DAY

    for step in range(1, steps):
        (
            current_positions,
            current_velocities
        ) = verlet_step(
            current_positions,
            current_velocities,
            masses,
            dt,
            softening
        )

        history[step] = current_positions

    return history


def centre_of_mass(positions, masses):
    return np.sum(
        positions * masses[:, None],
        axis=0
    ) / np.sum(masses)


def diagnostics(history, masses):
    initial = history[0]
    final = history[-1]

    initial_com = centre_of_mass(
        initial,
        masses
    )

    final_com = centre_of_mass(
        final,
        masses
    )

    displacement = np.linalg.norm(
        final_com - initial_com
    )

    print("\nDiagnostics")
    print("-" * 45)

    print(
        f"Initial COM: "
        f"{initial_com / AU}"
    )

    print(
        f"Final COM:   "
        f"{final_com / AU}"
    )

    print(
        f"COM displacement: "
        f"{displacement / AU:.6e} AU"
    )


def sun_earth():
    masses = np.array([
        1.0 * MSUN,
        3.003e-6 * MSUN
    ])

    positions = np.array([
        [0.0, 0.0, 0.0],
        [AU, 0.0, 0.0]
    ])

    earth_speed = np.sqrt(
        G * masses[0] / AU
    )

    velocities = np.array([
        [0.0, 0.0, 0.0],
        [0.0, earth_speed, 0.0]
    ])

    total_mass = np.sum(masses)

    velocities[0] = -(
        masses[1] / total_mass
    ) * velocities[1]

    return (
        positions,
        velocities,
        masses,
        ["Sun", "Earth"]
    )


def sun_earth_mars():
    masses = np.array([
        1.0 * MSUN,
        3.003e-6 * MSUN,
        3.227e-7 * MSUN
    ])

    positions = np.array([
        [0.0, 0.0, 0.0],
        [AU, 0.0, 0.0],
        [1.524 * AU, 0.0, 0.0]
    ])

    earth_speed = np.sqrt(
        G * masses[0] / AU
    )

    mars_speed = np.sqrt(
        G * masses[0] / (1.524 * AU)
    )

    velocities = np.array([
        [0.0, 0.0, 0.0],
        [0.0, earth_speed, 0.0],
        [0.0, mars_speed, 0.0]
    ])

    total_mass = np.sum(masses)

    velocities[0] = -(
        masses[1] * velocities[1]
        + masses[2] * velocities[2]
    ) / total_mass

    return (
        positions,
        velocities,
        masses,
        ["Sun", "Earth", "Mars"]
    )


def binary_stars():
    mass1 = 1.0 * MSUN
    mass2 = 0.8 * MSUN

    separation = 1.0 * AU

    r1 = (
        separation * mass2
        / (mass1 + mass2)
    )

    r2 = (
        separation * mass1
        / (mass1 + mass2)
    )

    positions = np.array([
        [-r1, 0.0, 0.0],
        [r2, 0.0, 0.0]
    ])

    angular_speed = np.sqrt(
        G * (mass1 + mass2)
        / separation**3
    )

    velocities = np.array([
        [0.0, -angular_speed * r1, 0.0],
        [0.0, angular_speed * r2, 0.0]
    ])

    masses = np.array([
        mass1,
        mass2
    ])

    return (
        positions,
        velocities,
        masses,
        ["Star A", "Star B"]
    )


def custom_system(n):
    masses = []
    positions = []
    velocities = []
    names = []

    print("\nEnter the parameters for each body.")
    print("Mass: kg")
    print("Radius: AU")
    print("Angles: degrees")
    print("Speed: km/s")

    for i in range(n):
        print(f"\nBody {i + 1}")

        mass = get_number(
            "Mass (kg): ",
            0
        )

        radius = get_number(
            "Radius (AU): ",
            0
        )

        angle = get_number(
            "Position angle (degrees): "
        )

        speed = get_number(
            "Speed (km/s): ",
            0
        )

        direction = get_number(
            "Velocity direction (degrees): "
        )

        masses.append(mass)

        positions.append(
            position(
                radius,
                angle
            )
        )

        velocities.append(
            velocity(
                speed,
                direction
            )
        )

        names.append(
            f"Body {i + 1}"
        )

    return (
        np.array(positions),
        np.array(velocities),
        np.array(masses),
        names
    )


def random_star_cluster(n=1000):
    rng = np.random.default_rng(42)

    masses = (
        rng.uniform(0.1, 2.0, n)
        * MSUN
    )

    positions = rng.normal(
        0,
        0.6 * AU,
        (n, 3)
    )

    velocities = rng.normal(
        0,
        5000,
        (n, 3)
    )

    total_mass = np.sum(masses)

    mean_velocity = np.sum(
        velocities * masses[:, None],
        axis=0
    ) / total_mass

    velocities -= mean_velocity

    com = centre_of_mass(
        positions,
        masses
    )

    positions -= com

    names = [
        f"Star {i + 1}"
        for i in range(n)
    ]

    return (
        positions,
        velocities,
        masses,
        names
    )


def body_style(name):
    if name == "Sun":
        return {
            "color": "gold",
            "marker": "o",
            "size": 260
        }

    if name == "Earth":
        return {
            "color": "royalblue",
            "marker": "o",
            "size": 70
        }

    if name == "Mars":
        return {
            "color": "red",
            "marker": "o",
            "size": 70
        }

    if name == "Star A":
        return {
            "color": "royalblue",
            "marker": "o",
            "size": 100
        }

    if name == "Star B":
        return {
            "color": "orange",
            "marker": "o",
            "size": 100
        }

    return {
        "color": "black",
        "marker": "o",
        "size": 60
    }


def create_legend(names):
    handles = []

    for name in names:
        style = body_style(name)

        handles.append(
            Line2D(
                [0],
                [0],
                marker=style["marker"],
                color="none",
                markerfacecolor=style["color"],
                markeredgecolor="black",
                markersize=8,
                label=name
            )
        )

    return handles


def plot_2d(history, names):
    fig, ax = plt.subplots(
        figsize=(10, 8)
    )

    for i, name in enumerate(names):
        style = body_style(name)

        x = history[:, i, 0] / AU
        y = history[:, i, 1] / AU

        ax.plot(
            x,
            y,
            linewidth=1.5
        )

        ax.scatter(
            x[0],
            y[0],
            s=style["size"],
            c=style["color"],
            marker=style["marker"],
            edgecolors="black",
            zorder=10
        )

    ax.set_xlabel(
        "x position (AU)",
        labelpad=8
    )

    ax.set_ylabel(
        "y position (AU)",
        labelpad=8
    )

    ax.set_title(
        "Computational N-Body Gravity Simulator - 2D"
    )

    ax.set_aspect("equal")
    ax.grid(True)

    ax.legend(
        handles=create_legend(names),
        loc="best"
    )

    fig.tight_layout()

    fig.savefig(
        "outputs/trajectory_2d.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


def plot_3d(history, names):
    fig = plt.figure(
        figsize=(10, 8)
    )

    ax = fig.add_subplot(
        111,
        projection="3d"
    )

    for i, name in enumerate(names):
        style = body_style(name)

        x = history[:, i, 0] / AU
        y = history[:, i, 1] / AU
        z = history[:, i, 2] / AU

        ax.plot(
            x,
            y,
            z,
            linewidth=1.5
        )

        ax.scatter(
            x[0],
            y[0],
            z[0],
            s=style["size"],
            c=style["color"],
            marker=style["marker"],
            edgecolors="black"
        )

    ax.set_xlabel(
        "x position (AU)",
        labelpad=10
    )

    ax.set_ylabel(
        "y position (AU)",
        labelpad=10
    )

    ax.set_zlabel(
        "z position (AU)",
        labelpad=10
    )

    ax.set_title(
        "Computational N-Body Gravity Simulator - 3D",
        pad=15
    )

    ax.legend(
        handles=create_legend(names),
        loc="best"
    )

    fig.tight_layout()

    fig.savefig(
        "outputs/trajectory_3d.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


def animate_2d(history, names):
    fig, ax = plt.subplots(
        figsize=(10, 8)
    )

    all_x = history[:, :, 0] / AU
    all_y = history[:, :, 1] / AU

    limit = max(
        np.max(np.abs(all_x)),
        np.max(np.abs(all_y))
    )

    if limit == 0:
        limit = 1

    limit *= 1.1

    ax.set_xlim(
        -limit,
        limit
    )

    ax.set_ylim(
        -limit,
        limit
    )

    ax.set_xlabel(
        "x position (AU)",
        labelpad=8
    )

    ax.set_ylabel(
        "y position (AU)",
        labelpad=8
    )

    ax.set_title(
        "Computational N-Body Gravity Simulator - Animation"
    )

    ax.set_aspect("equal")
    ax.grid(True)

    lines = []
    points = []

    for i, name in enumerate(names):
        style = body_style(name)

        line, = ax.plot(
            [],
            [],
            linewidth=1.5
        )

        point, = ax.plot(
            [],
            [],
            marker=style["marker"],
            linestyle="None",
            markersize=np.sqrt(
                style["size"]
            ),
            markerfacecolor=style["color"],
            markeredgecolor="black"
        )

        lines.append(line)
        points.append(point)

    ax.legend(
        handles=create_legend(names),
        loc="best"
    )

    def update(frame):
        for i in range(len(names)):
            x = (
                history[
                    :frame + 1,
                    i,
                    0
                ] / AU
            )

            y = (
                history[
                    :frame + 1,
                    i,
                    1
                ] / AU
            )

            lines[i].set_data(
                x,
                y
            )

            points[i].set_data(
                [x[-1]],
                [y[-1]]
            )

        return lines + points

    animation = FuncAnimation(
        fig,
        update,
        frames=len(history),
        interval=40,
        repeat=False
    )

    animation.save(
        "outputs/n_body_animation.gif",
        writer=PillowWriter(fps=25)
    )

    plt.show()


def plot_cluster_2d(history):
    final_positions = (
        history[-1] / AU
    )

    fig, ax = plt.subplots(
        figsize=(10, 8)
    )

    ax.set_facecolor("black")

    ax.scatter(
        final_positions[:, 0],
        final_positions[:, 1],
        s=3,
        c="red"
    )

    ax.scatter(
        [0],
        [0],
        s=180,
        c="yellow",
        edgecolors="white",
        zorder=10,
        label="Cluster centre"
    )

    ax.set_xlabel(
        "x position (AU)",
        color="white",
        labelpad=8
    )

    ax.set_ylabel(
        "y position (AU)",
        color="white",
        labelpad=8
    )

    ax.set_title(
        "1000-Star N-Body Cluster - 2D",
        color="white"
    )

    ax.tick_params(
        colors="white"
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        "outputs/star_cluster_1000_2d.png",
        dpi=300,
        facecolor="black",
        bbox_inches="tight"
    )

    plt.show()


def plot_cluster_3d(history):
    final_positions = (
        history[-1] / AU
    )

    fig = plt.figure(
        figsize=(10, 8)
    )

    ax = fig.add_subplot(
        111,
        projection="3d"
    )

    ax.scatter(
        final_positions[:, 0],
        final_positions[:, 1],
        final_positions[:, 2],
        s=3,
        c="red"
    )

    ax.set_xlabel(
        "x position (AU)",
        labelpad=10
    )

    ax.set_ylabel(
        "y position (AU)",
        labelpad=10
    )

    ax.set_zlabel(
        "z position (AU)",
        labelpad=10
    )

    ax.set_title(
        "1000-Star N-Body Cluster - 3D",
        pad=15
    )

    fig.tight_layout()

    fig.savefig(
        "outputs/star_cluster_1000_3d.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


def animate_cluster(history):
    fig, ax = plt.subplots(
        figsize=(10, 8)
    )

    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    positions = history / AU

    limit = max(
        np.max(np.abs(positions[:, :, 0])),
        np.max(np.abs(positions[:, :, 1]))
    )

    if limit == 0:
        limit = 1

    limit *= 1.1

    ax.set_xlim(
        -limit,
        limit
    )

    ax.set_ylim(
        -limit,
        limit
    )

    ax.set_xlabel(
        "x position (AU)",
        color="white",
        labelpad=8
    )

    ax.set_ylabel(
        "y position (AU)",
        color="white",
        labelpad=8
    )

    ax.set_title(
        "1000-Star N-Body Cluster",
        color="white"
    )

    ax.tick_params(
        colors="white"
    )

    scatter = ax.scatter(
        [],
        [],
        s=4,
        c="red"
    )

    def update(frame):
        xy = positions[
            frame,
            :,
            :2
        ]

        scatter.set_offsets(xy)

        return scatter,

    animation = FuncAnimation(
        fig,
        update,
        frames=len(history),
        interval=40,
        repeat=False
    )

    animation.save(
        "outputs/star_cluster_1000.gif",
        writer=PillowWriter(fps=25)
    )

    plt.show()


def comparison_data():
    systems = [
        sun_earth,
        sun_earth_mars,
        binary_stars
    ]

    results = []

    for system in systems:
        positions, velocities, masses, names = system()

        history = simulate(
            positions,
            velocities,
            masses,
            365.25,
            2.0
        )

        results.append(
            (
                history,
                names
            )
        )

    return results


def comparison_2d():
    results = comparison_data()

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18, 7)
    )

    titles = [
        "Sun-Earth",
        "Sun-Earth-Mars",
        "Binary Stars"
    ]

    for ax, result, title in zip(
        axes,
        results,
        titles
    ):
        history, names = result

        for i, name in enumerate(names):
            style = body_style(name)

            x = history[:, i, 0] / AU
            y = history[:, i, 1] / AU

            ax.plot(
                x,
                y,
                linewidth=1.5
            )

            ax.scatter(
                x[0],
                y[0],
                s=style["size"],
                c=style["color"],
                marker=style["marker"],
                edgecolors="black",
                zorder=10
            )

        ax.set_title(
            title,
            pad=12
        )

        ax.set_xlabel(
            "x position (AU)",
            labelpad=8
        )

        ax.set_ylabel(
            "y position (AU)",
            labelpad=8
        )

        ax.grid(True)
        ax.set_aspect("equal")

        ax.legend(
            handles=create_legend(names),
            loc="best"
        )

    fig.suptitle(
        "Three-System N-Body Gravity Comparison",
        fontsize=16,
        y=0.98
    )

    fig.subplots_adjust(
        left=0.06,
        right=0.98,
        bottom=0.12,
        top=0.86,
        wspace=0.28
    )

    fig.savefig(
        "outputs/comparison_2d.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


def comparison_3d():
    results = comparison_data()

    fig = plt.figure(
        figsize=(19, 7)
    )

    titles = [
        "Sun-Earth",
        "Sun-Earth-Mars",
        "Binary Stars"
    ]

    for index, (result, title) in enumerate(
        zip(results, titles),
        start=1
    ):
        history, names = result

        ax = fig.add_subplot(
            1,
            3,
            index,
            projection="3d"
        )

        for i, name in enumerate(names):
            style = body_style(name)

            x = history[:, i, 0] / AU
            y = history[:, i, 1] / AU
            z = history[:, i, 2] / AU

            ax.plot(
                x,
                y,
                z,
                linewidth=1.5
            )

            ax.scatter(
                x[0],
                y[0],
                z[0],
                s=style["size"],
                c=style["color"],
                marker=style["marker"],
                edgecolors="black"
            )

        ax.set_title(
            title,
            pad=15
        )

        ax.set_xlabel(
            "x position (AU)",
            labelpad=10
        )

        ax.set_ylabel(
            "y position (AU)",
            labelpad=10
        )

        ax.set_zlabel(
            "z position (AU)",
            labelpad=10
        )

        ax.legend(
            handles=create_legend(names),
            loc="best"
        )

    fig.suptitle(
        "Three-System N-Body Gravity Comparison",
        fontsize=16,
        y=0.98
    )

    fig.subplots_adjust(
        left=0.02,
        right=0.98,
        bottom=0.08,
        top=0.84,
        wspace=0.15
    )

    fig.savefig(
        "outputs/comparison_3d.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


def comparison_animation():
    results = comparison_data()

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18, 7)
    )

    titles = [
        "Sun-Earth",
        "Sun-Earth-Mars",
        "Binary Stars"
    ]

    lines = []
    points = []

    for ax, result, title in zip(
        axes,
        results,
        titles
    ):
        history, names = result

        positions = history[:, :, :2] / AU

        limit = max(
            np.max(np.abs(positions[:, :, 0])),
            np.max(np.abs(positions[:, :, 1]))
        )

        if limit == 0:
            limit = 1

        limit *= 1.1

        ax.set_xlim(
            -limit,
            limit
        )

        ax.set_ylim(
            -limit,
            limit
        )

        ax.set_title(
            title,
            pad=12
        )

        ax.set_xlabel(
            "x position (AU)",
            labelpad=8
        )

        ax.set_ylabel(
            "y position (AU)",
            labelpad=8
        )

        ax.grid(True)
        ax.set_aspect("equal")

        system_lines = []
        system_points = []

        for i, name in enumerate(names):
            style = body_style(name)

            line, = ax.plot(
                [],
                [],
                linewidth=1.5
            )

            point, = ax.plot(
                [],
                [],
                marker=style["marker"],
                linestyle="None",
                markersize=np.sqrt(
                    style["size"]
                ),
                markerfacecolor=style["color"],
                markeredgecolor="black"
            )

            system_lines.append(line)
            system_points.append(point)

        lines.append(system_lines)
        points.append(system_points)

        ax.legend(
            handles=create_legend(names),
            loc="best"
        )

    fig.suptitle(
        "Three-System N-Body Gravity Comparison",
        fontsize=16,
        y=0.98
    )

    fig.subplots_adjust(
        left=0.06,
        right=0.98,
        bottom=0.12,
        top=0.86,
        wspace=0.28
    )

    frames = min(
        len(results[0][0]),
        len(results[1][0]),
        len(results[2][0])
    )

    def update(frame):
        artists = []

        for system_index, result in enumerate(
            results
        ):
            history, names = result

            for i in range(len(names)):
                x = (
                    history[
                        :frame + 1,
                        i,
                        0
                    ] / AU
                )

                y = (
                    history[
                        :frame + 1,
                        i,
                        1
                    ] / AU
                )

                lines[
                    system_index
                ][i].set_data(
                    x,
                    y
                )

                points[
                    system_index
                ][i].set_data(
                    [x[-1]],
                    [y[-1]]
                )

                artists.append(
                    lines[system_index][i]
                )

                artists.append(
                    points[system_index][i]
                )

        return artists

    animation = FuncAnimation(
        fig,
        update,
        frames=frames,
        interval=50,
        repeat=False
    )

    animation.save(
        "outputs/three_system_comparison.gif",
        writer=PillowWriter(fps=25)
    )

    plt.show()


def run_standard_system(system_function):
    print("\nSelect visualization:")
    print("1. 2D")
    print("2. 3D")
    print("3. Animated 2D")

    visualization = input(
        "Enter choice: "
    ).strip()

    duration = get_number(
        "Simulation duration (days): ",
        0
    )

    dt = get_number(
        "Time step (days): ",
        0
    )

    positions, velocities, masses, names = system_function()

    history = simulate(
        positions,
        velocities,
        masses,
        duration,
        dt
    )

    diagnostics(
        history,
        masses
    )

    if visualization == "1":
        plot_2d(
            history,
            names
        )

    elif visualization == "2":
        plot_3d(
            history,
            names
        )

    elif visualization == "3":
        animate_2d(
            history,
            names
        )

    else:
        print("Invalid visualization choice.")


def main():
    print("=" * 60)
    print("COMPUTATIONAL N-BODY GRAVITY SIMULATOR")
    print("=" * 60)

    print("\nSelect a system:")
    print("1. Sun - Earth")
    print("2. Sun - Earth - Mars")
    print("3. Binary Stars")
    print("4. Custom N-body System")
    print("5. Automatic 1000-Star N-body Cluster")
    print("6. Three-System Comparison")

    choice = input(
        "\nEnter choice: "
    ).strip()

    if choice == "1":
        run_standard_system(
            sun_earth
        )

    elif choice == "2":
        run_standard_system(
            sun_earth_mars
        )

    elif choice == "3":
        run_standard_system(
            binary_stars
        )

    elif choice == "4":
        n = int(
            get_number(
                "Number of bodies: ",
                2
            )
        )

        positions, velocities, masses, names = custom_system(
            n
        )

        print("\nSelect visualization:")
        print("1. 2D")
        print("2. 3D")
        print("3. Animated 2D")

        visualization = input(
            "Enter choice: "
        ).strip()

        duration = get_number(
            "Simulation duration (days): ",
            0
        )

        dt = get_number(
            "Time step (days): ",
            0
        )

        history = simulate(
            positions,
            velocities,
            masses,
            duration,
            dt
        )

        diagnostics(
            history,
            masses
        )

        if visualization == "1":
            plot_2d(
                history,
                names
            )

        elif visualization == "2":
            plot_3d(
                history,
                names
            )

        elif visualization == "3":
            animate_2d(
                history,
                names
            )

        else:
            print("Invalid visualization choice.")

    elif choice == "5":
        print("\nAutomatic 1000-Star N-body Cluster")

        duration = get_number(
            "Simulation duration (days): ",
            0
        )

        dt = get_number(
            "Time step (days): ",
            0
        )

        softening_au = get_number(
            "Softening length (AU): ",
            0
        )

        positions, velocities, masses, names = random_star_cluster(
            1000
        )

        history = simulate(
            positions,
            velocities,
            masses,
            duration,
            dt,
            softening_au * AU
        )

        print("\nSelect visualization:")
        print("1. 2D Black/Red")
        print("2. 3D")
        print("3. Animated 2D Black/Red")

        visualization = input(
            "Enter choice: "
        ).strip()

        if visualization == "1":
            plot_cluster_2d(
                history
            )

        elif visualization == "2":
            plot_cluster_3d(
                history
            )

        elif visualization == "3":
            animate_cluster(
                history
            )

        else:
            print("Invalid visualization choice.")

    elif choice == "6":
        print("\nThree-System Comparison")

        print("1. 2D comparison")
        print("2. 3D comparison")
        print("3. Animated comparison")

        visualization = input(
            "Enter choice: "
        ).strip()

        if visualization == "1":
            comparison_2d()

        elif visualization == "2":
            comparison_3d()

        elif visualization == "3":
            comparison_animation()

        else:
            print("Invalid visualization choice.")

    else:
        print("Invalid system choice.")


if __name__ == "__main__":
    main()