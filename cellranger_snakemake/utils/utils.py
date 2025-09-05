# Function to apply suffix if specified
def get_directories(config):
    """Get directory paths with optional suffix applied"""
    dirs = {}
    base_dirs = config.get("directories")
    suffix = config.get("directories_suffix", "none")
    
    for key, value in base_dirs.items():
        if suffix != "none" and suffix != "":
            # Apply suffix to directory path
            dirs[key] = f"{value}_{suffix}"
        else:
            # Use directory path as-is
            dirs[key] = value
    
    return dirs

def has_underscore(s: str) -> bool:
    """
    Return True if the string contains at least one underscore, else False.
    """
    return "_" in s