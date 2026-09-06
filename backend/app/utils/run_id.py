import time
import random
import string
from datetime import datetime

def new_run_id() -> str:
    """
    Generate a collision-safe run ID suitable for filesystem directories.
    Format: YYYYMMDD-HHMM-{6char}
    Example: 20250906-1015-A3X9KL
    """
    # Get current timestamp
    now = datetime.now()
    date_part = now.strftime("%Y%m%d-%H%M")
    
    # Generate 6-character random suffix (alphanumeric, uppercase)
    chars = string.ascii_uppercase + string.digits
    random_part = ''.join(random.choices(chars, k=6))
    
    return f"{date_part}-{random_part}"

def sanitize_run_id(run_id: str) -> str:
    """
    Sanitize run_id to be filesystem-safe.
    Removes/replaces problematic characters.
    """
    # Replace problematic characters with underscores
    problematic = ['/', '\\', ':', '*', '?', '"', '<', '>', '|', ' ']
    sanitized = run_id
    for char in problematic:
        sanitized = sanitized.replace(char, '_')
    
    # Limit length to reasonable filesystem limits
    if len(sanitized) > 100:
        sanitized = sanitized[:100]
    
    return sanitized