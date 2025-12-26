import os

EXCLUDE_DIRS = {".venv"}

def print_tree(path, prefix=""):
    """Рекурсивно виводить дерево папок і файлів,
    ігноруючи каталоги з EXCLUDE_DIRS.
    """
    try:
        items = sorted(os.listdir(path))
    except PermissionError:
        return

    # фільтруємо виключені папки
    items = [
        item for item in items
        if item not in EXCLUDE_DIRS
    ]

    for index, name in enumerate(items):
        full_path = os.path.join(path, name)
        connector = "└── " if index == len(items) - 1 else "├── "
        print(prefix + connector + name)

        if os.path.isdir(full_path):
            new_prefix = prefix + ("    " if index == len(items) - 1 else "│   ")
            print_tree(full_path, new_prefix)


# 🔹 Корінь проєкту (поточна директорія)
project_root = os.getcwd()

print(f"Дерево проєкту: {project_root}\n")
print_tree(project_root)