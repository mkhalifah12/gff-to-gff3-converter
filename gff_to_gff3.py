from Bio import SeqIO

def convert_gff_to_gff3(input_gff, output_gff3):
    """
    Converts a GFF file to GFF3 format using Biopython.
    
    Args:
        input_gff (str): Path to the input GFF file.
        output_gff3 (str): Path to the output GFF3 file.
    """
    with open(input_gff, "r") as input_handle, open(output_gff3, "w") as output_handle:
        # Write the GFF3 header
        output_handle.write("##gff-version 3\n")
        
        # Track parent-child relationships
        feature_dict = {}  # Stores parent features and their children
        
        for line in input_handle:
            # Skip comment lines
            if line.startswith("#"):
                continue
            
            # Split the GFF line into columns
            columns = line.strip().split("\t")
            if len(columns) != 9:
                continue  # Skip malformed lines
            
            # Extract attributes
            attributes = columns[8]
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attr_dict[key] = value
                elif " " in attr:
                    key, value = attr.split(" ", 1)
                    attr_dict[key] = value.strip('"')
            
            # Handle ID and Parent attributes
            if "ID" not in attr_dict:
                # Generate a unique ID if missing
                feature_type = columns[2]
                attr_dict["ID"] = f"{feature_type}_{columns[3]}_{columns[4]}"
            
            if "Parent" not in attr_dict and columns[2] != "gene":
                # Assign Parent attribute for child features (e.g., exons, CDS)
                parent_id = f"gene_{columns[3]}_{columns[4]}"
                attr_dict["Parent"] = parent_id
            
            # Rebuild the attributes column
            new_attributes = ";".join([f"{key}={value}" for key, value in attr_dict.items()])
            columns[8] = new_attributes
            
            # Write the modified line to the output file
            output_handle.write("\t".join(columns) + "\n")

if __name__ == "__main__":
    # Input and output file paths
    input_gff = "input.gff"  # Replace with your input GFF file path
    output_gff3 = "output.gff3"  # Replace with your output GFF3 file path
    
    # Run the conversion
    convert_gff_to_gff3(input_gff, output_gff3)
    print(f"Conversion complete! Output saved to {output_gff3}")