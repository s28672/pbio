import random
import os
import re  # Import for input validation


def generate_dna_sequence(length):
    """Generate a random DNA sequence of specified length"""
    # Define the possible nucleotides in DNA
    nucleotides = ['A', 'C', 'G', 'T']
    # Generate a random sequence of specified length
    sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    return sequence


def insert_name_into_sequence(sequence, name):
    """Insert the name at a random position in the sequence"""
    # Choose a random position to insert the name
    position = random.randint(0, len(sequence))
    # Insert the name at the chosen position
    modified_sequence = sequence[:position] + name + sequence[position:]
    return modified_sequence


def calculate_statistics(sequence, name):
    """Calculate statistics of nucleotides in the sequence (excluding the name)"""
    # Remove the name from the sequence for statistics calculation
    pure_sequence = sequence.replace(name, "")

    # Calculate the count of each nucleotide
    total_length = len(pure_sequence)
    a_count = pure_sequence.count('A')
    c_count = pure_sequence.count('C')
    g_count = pure_sequence.count('G')
    t_count = pure_sequence.count('T')

    # Calculate percentages
    a_percent = (a_count / total_length) * 100
    c_percent = (c_count / total_length) * 100
    g_percent = (g_count / total_length) * 100
    t_percent = (t_count / total_length) * 100

    # ORIGINAL:
    # # Calculate CG ratio (% of C and G nucleotides to A and T)
    # cg_percent = ((c_count + g_count) / total_length) * 100
    # MODIFIED (fixed the ratio calculation to be C+G to A+T):
    cg_ratio = (c_count + g_count) / (a_count + t_count) if (a_count + t_count) > 0 else 0
    cg_percent = ((c_count + g_count) / total_length) * 100

    return {
        'A': a_percent,
        'C': c_percent,
        'G': g_percent,
        'T': t_percent,
        'CG_percent': cg_percent,
        'CG_ratio': cg_ratio
    }


def save_to_fasta(sequence_id, description, sequence):
    """Save the sequence to a FASTA format file"""
    # ORIGINAL:
    # # Create the filename based on the sequence ID
    # filename = f"{sequence_id}.fasta"
    # MODIFIED (sanitize filename to prevent invalid characters):
    # Sanitize the sequence_id to create a valid filename
    safe_id = re.sub(r'[^\w\-\.]', '_', sequence_id)
    filename = f"{safe_id}.fasta"

    # Create the FASTA header
    header = f">{sequence_id} {description}"

    # Write to the file
    with open(filename, 'w') as file:
        file.write(header + '\n')
        # ORIGINAL:
        # file.write(sequence + '\n')
        # MODIFIED (format sequence with line breaks for better readability):
        # Format sequence with line breaks (standard FASTA format uses 60-80 chars per line)
        for i in range(0, len(sequence), 60):
            file.write(sequence[i:i + 60] + '\n')

    return filename


def main():
    """Main function to run the DNA sequence generator"""
    print("DNA Sequence Generator in FASTA Format")
    print("--------------------------------------")

    # Get user input for sequence length
    while True:
        try:
            length = int(input("Enter the sequence length: "))
            if length <= 0:
                print("Please enter a positive number.")
                continue
            break
        except ValueError:
            print("Please enter a valid number.")

    # Get sequence ID and description
    while True:
        sequence_id = input("Enter the sequence ID: ")
        if sequence_id.strip():  # Check if ID is not empty
            break
        print("Sequence ID cannot be empty. Please try again.")

    description = input("Provide a description of the sequence: ")

    # Get user's name
    while True:
        name = input("Enter your name: ")
        if name.strip():  # Check if name is not empty
            break
        print("Name cannot be empty. Please try again.")

    # Generate the sequence
    sequence = generate_dna_sequence(length)

    # Insert name into the sequence
    modified_sequence = insert_name_into_sequence(sequence, name)

    # Save to FASTA file
    filename = save_to_fasta(sequence_id, description, modified_sequence)

    # Calculate statistics
    stats = calculate_statistics(modified_sequence, name)

    # Display results
    print(f"\nThe sequence was saved to the file {filename}")
    print("Sequence statistics:")
    print(f"A: {stats['A']:.1f}%")
    print(f"C: {stats['C']:.1f}%")
    print(f"G: {stats['G']:.1f}%")
    print(f"T: {stats['T']:.1f}%")
    # ORIGINAL:
    # print(f"%CG: {stats['CG']:.1f}")
    # MODIFIED (added both CG percentage and ratio for clarity):
    print(f"%CG: {stats['CG_percent']:.1f}")
    print(f"CG/AT ratio: {stats['CG_ratio']:.2f}")


if __name__ == "__main__":
    main()
