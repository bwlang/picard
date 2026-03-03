package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Utility class for reading interval files in BED or Picard interval_list format.
 *
 * <p>This class provides format-detection, generic loading (auto-detecting BED vs interval_list),
 * and a dedicated BED parser.  It is designed to be used by any Picard tool that needs to
 * accept interval files — callers get transparent support for both formats without duplicating
 * parsing logic.
 *
 * <p>BED coordinate conventions: BED uses 0-based half-open coordinates.  This parser converts
 * them to the 1-based closed intervals used throughout Picard / htsjdk.
 */
public final class IntervalFileReader {

    /** Recognized interval file formats. */
    public enum IntervalFileFormat { INTERVAL_LIST, BED, UNKNOWN }

    /** Result of format detection: the detected format and, when format is UNKNOWN, the first data line. */
    public record FormatDetectionResult(IntervalFileFormat format, String firstLine) {}

    private IntervalFileReader() {} // utility class — no instances

    /**
     * Detects whether a reader contains interval_list or BED content by inspecting the first
     * significant (non-empty, non-comment) line.  The reader's position is NOT reset after this
     * call — callers must {@link BufferedReader#mark} and {@link BufferedReader#reset} as needed.
     *
     * <ul>
     *   <li>Lines starting with {@code @} indicate interval_list (SAM-style header).</li>
     *   <li>Lines with ≥3 tab-separated fields indicate BED.</li>
     *   <li>Comment lines starting with {@code #} are skipped.</li>
     * </ul>
     */
    public static FormatDetectionResult detectIntervalFormat(final BufferedReader reader) throws IOException {
        String line;
        while ((line = reader.readLine()) != null) {
            final String trimmed = line.trim();
            if (trimmed.isEmpty() || trimmed.startsWith("#")) {
                continue;
            }
            if (trimmed.startsWith("@")) {
                return new FormatDetectionResult(IntervalFileFormat.INTERVAL_LIST, null);
            } else if (trimmed.split("\t").length >= 3) {
                return new FormatDetectionResult(IntervalFileFormat.BED, null);
            } else {
                return new FormatDetectionResult(IntervalFileFormat.UNKNOWN, trimmed);
            }
        }
        return new FormatDetectionResult(IntervalFileFormat.UNKNOWN, null);
    }

    /**
     * Loads an {@link IntervalList} from {@code reader}, auto-detecting whether the content is
     * BED or interval_list format.  The reader must support {@link BufferedReader#mark} /
     * {@link BufferedReader#reset} (a plain {@link BufferedReader} always does — this works for
     * regular files, pipes, FIFOs, and process-substitution streams alike).
     *
     * <p>Returns a sorted, uniqued {@link IntervalList}.
     *
     * @param reader     source; must support mark/reset
     * @param dictionary sequence dictionary used when parsing BED-format content
     * @throws IOException     on read error
     * @throws PicardException on unrecognized format
     */
    public static IntervalList loadIntervals(final BufferedReader reader,
                                             final SAMSequenceDictionary dictionary) throws IOException {
        // 8 KB is enough to cover any realistic BED/interval_list preamble.
        // BufferedReader.mark() is backed by an in-memory char array, so this also works
        // correctly for non-seekable sources such as pipes and FIFOs.
        reader.mark(8 * 1024);
        final FormatDetectionResult detected = detectIntervalFormat(reader);
        reader.reset();
        return switch (detected.format()) {
            case INTERVAL_LIST -> IntervalList.fromReader(reader).uniqued();
            case BED -> {
                final SAMFileHeader header = new SAMFileHeader();
                header.setSequenceDictionary(dictionary);
                yield fromBed(reader, header, false, false, null).uniqued();
            }
            case UNKNOWN -> throw new PicardException(
                    "Unrecognized interval file format. Expected interval_list (lines starting with @) " +
                    "or BED (≥3 tab-separated fields). First data line: " + detected.firstLine());
        };
    }

    /**
     * Convenience overload of {@link #loadIntervals(BufferedReader, SAMSequenceDictionary)} that
     * opens {@code file} itself.  {@link IOException} is wrapped in a {@link PicardException}.
     *
     * @param file       BED or interval_list file to read
     * @param dictionary sequence dictionary used when parsing BED-format content
     */
    public static IntervalList loadIntervals(final File file, final SAMSequenceDictionary dictionary) {
        try (final BufferedReader reader = new BufferedReader(new FileReader(file))) {
            return loadIntervals(reader, dictionary);
        } catch (final IOException e) {
            throw new PicardException("Error reading intervals from " + file, e);
        }
    }

    /**
     * Parses BED records from a {@link BufferedReader} into an {@link IntervalList}.
     * Comment lines starting with {@code #} are skipped.
     * Coordinates are converted from 0-based half-open BED to 1-based closed intervals.
     * <p>
     * This method works with any reader, including streams backed by pipes or FIFOs.
     *
     * @param reader                  source of BED records
     * @param header                  SAMFileHeader whose sequence dictionary is used for validation
     * @param dropMissingContigs      if true, records on contigs absent from the dictionary are silently
     *                                dropped; if false, such records throw a {@link PicardException}
     * @param keepLengthZeroIntervals if true, length-zero BED intervals (chromStart == chromEnd) are
     *                                included; if false, they are silently dropped
     * @param log                     used for informational and warning messages; may be {@code null}
     *                                to suppress all logging
     * @return an unsorted, non-uniqued {@link IntervalList}
     * @throws IOException on read error
     */
    public static IntervalList fromBed(final BufferedReader reader,
                                       final SAMFileHeader header,
                                       final boolean dropMissingContigs,
                                       final boolean keepLengthZeroIntervals,
                                       final Log log) throws IOException {
        final IntervalList intervalList = new IntervalList(header);
        int missingIntervals = 0;
        int missingRegion = 0;
        int lengthZeroIntervals = 0;

        String line;
        while ((line = reader.readLine()) != null) {
            final String trimmed = line.trim();
            if (trimmed.isEmpty() || trimmed.startsWith("#")) {
                continue;
            }

            final String[] fields = trimmed.split("\t");
            if (fields.length < 3) {
                throw new PicardException("Invalid BED line (fewer than 3 tab-separated fields): " + line);
            }

            final String sequenceName = fields[0];
            final int start = Integer.parseInt(fields[1]) + 1; // BED start is 0-based; convert to 1-based
            final int end   = Integer.parseInt(fields[2]);     // BED end is 0-based exclusive == 1-based inclusive
            final String name = (fields.length > 3 && !fields[3].isEmpty()) ? fields[3] : null;
            final boolean isNegativeStrand = fields.length > 5 && "-".equals(fields[5]);

            final SAMSequenceRecord sequenceRecord = header.getSequenceDictionary().getSequence(sequenceName);

            if (null == sequenceRecord) {
                if (dropMissingContigs) {
                    if (log != null) log.info(String.format("Dropping interval with missing contig: %s:%d-%d", sequenceName, start, end));
                    missingIntervals++;
                    missingRegion += end - start + 1;
                    continue;
                }
                throw new PicardException(String.format("Sequence '%s' was not found in the sequence dictionary", sequenceName));
            } else if (start < 1) {
                throw new PicardException(String.format("Start on sequence '%s' was less than one: %d", sequenceName, start));
            } else if (sequenceRecord.getSequenceLength() < start) {
                throw new PicardException(String.format("Start on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), start));
            } else if ((end == 0 && start != 1) || end < 0) {
                throw new PicardException(String.format("End on sequence '%s' was less than one: %d", sequenceName, end));
            } else if (sequenceRecord.getSequenceLength() < end) {
                throw new PicardException(String.format("End on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), end));
            } else if (end < start - 1) {
                throw new PicardException(String.format("On sequence '%s', end < start-1: %d <= %d", sequenceName, end, start));
            }

            if ((start == end + 1) && !keepLengthZeroIntervals) {
                if (log != null) log.info(String.format("Skipping writing length zero interval at %s:%d-%d.", sequenceName, start, end));
                lengthZeroIntervals++;
                continue;
            }
            if (start == end + 1) {
                lengthZeroIntervals++;
            }

            intervalList.add(new Interval(sequenceName, start, end, isNegativeStrand, name));
        }

        if (log != null) {
            if (dropMissingContigs) {
                if (missingRegion == 0) {
                    log.info("There were no missing regions.");
                } else {
                    log.warn(String.format("There were %d missing regions with a total of %d bases", missingIntervals, missingRegion));
                }
            }
            if (!keepLengthZeroIntervals) {
                if (lengthZeroIntervals == 0) {
                    log.info("No input regions had length zero, so none were skipped.");
                } else {
                    log.info(String.format("Skipped writing a total of %d entries with length zero in the input file.", lengthZeroIntervals));
                }
            } else {
                if (lengthZeroIntervals > 0) {
                    log.warn(String.format("Input file had %d entries with length zero. Run with the KEEP_LENGTH_ZERO_INTERVALS flag set to false to remove these.", lengthZeroIntervals));
                }
            }
        }

        return intervalList;
    }
}
