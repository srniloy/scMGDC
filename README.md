# Installation Instructions

## Prerequisites
Ensure you have the following installed on your system:
- Python (>=3.8)
- `pip` (Python package manager)

## Setup Virtual Environment

1. **Clone the Repository**
   ```sh
   git clone <repository-url>
   cd <project-directory>
   ```

2. **Create a Virtual Environment**
   ```sh
   python -m venv venv
   ```

3. **Activate the Virtual Environment**
   - On macOS/Linux:
     ```sh
     source venv/bin/activate
     ```
   - On Windows:
     ```sh
     venv\Scripts\activate
     ```

4. **Install Dependencies**
   ```sh
   pip install -r requirements.txt
   ```

## Running the Project

1. **Ensure the Virtual Environment is Activated**
   ```sh
   source venv/bin/activate  # macOS/Linux
   venv\Scripts\activate    # Windows
   ```

2. **Run the Application**
   ```sh
   python scMGDC.py  # Replace with the actual entry point file
   ```

## Deactivating the Virtual Environment
To exit the virtual environment, run:
```sh
deactivate
```

## Additional Notes
- If `pip` is outdated, update it before installing dependencies:
  ```sh
  pip install --upgrade pip
  ```
- If facing issues, try reinstalling dependencies:
  ```sh
  pip install --no-cache-dir -r requirements.txt
  
